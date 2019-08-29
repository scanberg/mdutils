#pragma once

#include <core/common.h>
#include <core/types.h>
#include <core/array_types.h>

struct Bitfield {
	typedef uint32 ElementType;

	ElementType* ptr = nullptr;
	int64 count = 0;

	constexpr int64 size() const { return count; }
	constexpr int64 size_in_bytes() const { return (count + 8 - 1) / 8; } // @NOTE: round up integer div.
	
	constexpr ElementType* data() { return ptr; }
	constexpr ElementType* beg() { return ptr; }
	constexpr ElementType* end() { return ptr + count; }

	constexpr const ElementType* data() const { return ptr; }
	constexpr const ElementType* beg() const { return ptr; }
	constexpr const ElementType* end() const { return ptr + count; }

	operator bool() const {
		return ptr != nullptr;
	}
};

namespace bitfield {

namespace detail {
constexpr auto block_bits = sizeof(Bitfield::ElementType) * 8;

inline Bitfield::ElementType block(int64 idx) {
	return (Bitfield::ElementType)(idx / block_bits);
}

inline Bitfield::ElementType bit_pattern(int64 idx) {
	return (Bitfield::ElementType)1 << (idx % block_bits);
}

inline int64 number_of_set_bits(uint32 i) {
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

inline int64 num_blocks(Bitfield field) {
	return (field.size_in_bytes() + sizeof(Bitfield::ElementType) - 1) / sizeof(Bitfield::ElementType); // @NOTE: round up integer div.
}
}

inline void free(Bitfield* field) {
    ASSERT(field);
    if (field->ptr) {
        ALIGNED_FREE(field->ptr);
        field->ptr = nullptr;
        field->count = 0;
    }
}

inline void init(Bitfield* field, int64 num_bits) {
    ASSERT(field);
    free(field);
	field->count = num_bits;
    field->ptr = (Bitfield::ElementType*)ALIGNED_MALLOC(field->size_in_bytes(), 16);   // @NOTE: Align to 16 byte to allow for aligned simd load/store
	memset(field->ptr, 0, field->size_in_bytes());
}

inline void init(Bitfield* field, Bitfield src) {
	ASSERT(field);
	free(field);
	if (!src) return;
	field->ptr = (Bitfield::ElementType*)ALIGNED_MALLOC(src.size_in_bytes(), 16);
	field->count = src.count;
	memcpy(field->ptr, src.ptr, src.size_in_bytes());
}

inline void copy(Bitfield dst, const Bitfield src) {
	ASSERT(dst.size() == src.size() && "Bitfield size did not match");
	memcpy(dst.ptr, src.ptr, dst.size_in_bytes());
}

inline void set_all(Bitfield field) {
    memset(field.ptr, 0xFF, field.size_in_bytes());
}

inline void clear_all(Bitfield field) {
    memset(field.ptr, 0, field.size_in_bytes());
}

inline void invert_all(Bitfield field) {
    // @TODO: Vectorize 
    for (int64 i = 0; i < detail::num_blocks(field); i++) {
        field.ptr[i] = ~field.ptr[i];
    }
}

inline int64 number_of_bits_set(const Bitfield field) {
	const uint32* ptr = (uint32*)field.data();
	const int64 stride = sizeof(uint32) * 8;
	const int64 size = field.size() / stride;
	const int64 rest = field.size() % stride;
	int64 count = 0;
	for (int64 i = 0; i < size; i++) {
		count += detail::number_of_set_bits(ptr[i]);
	}
	if (rest != 0) {
		const auto bits = ptr[size] & (detail::bit_pattern(rest) - 1);
		count += detail::number_of_set_bits(bits);
	}

	return count;
}

template <typename Int>
inline void set_range(Bitfield field, Range<Int> range) {
	const auto beg_blk = detail::block(range.beg);
	const auto end_blk = detail::block(range.end);

	if (beg_blk == end_blk) {
		// All bits reside within the same Block
		const auto bits = (detail::bit_pattern(range.beg) - 1) ^ (detail::bit_pattern(range.end) - 1);
		field.ptr[beg_blk] |= bits;
		return;
	}

	field.ptr[beg_blk] |= (~(detail::bit_pattern(range.beg) - 1));
	field.ptr[end_blk] |= (detail::bit_pattern(range.end) - 1);

	// Set any bits within the inner range of blocks: beg_blk, [inner range], end_blk
	const int64 size = end_blk - beg_blk - 1;
	if (size > 0) {
		memset(field.ptr + beg_blk + 1, 0xFF, size * sizeof(Bitfield::ElementType));
	}
}

inline bool any_bit_set(const Bitfield field) {
	const auto size = field.size_in_bytes();
	const uint8* buf = (const uint8*)field.data();
	return !(buf[0] == 0 && !memcmp(buf, buf + 1, size - 1));
}

template <typename Int>
inline bool any_bit_set_in_range(const Bitfield field, Range<Int> range) {
	const auto beg_blk = detail::block(range.beg);
	const auto end_blk = detail::block(range.end);

	if (beg_blk == end_blk) {
		// All bits reside within the same Block
		const auto bit_mask = (detail::bit_pattern(range.beg) - 1) ^ (detail::bit_pattern(range.end) - 1);
		return (field.ptr[beg_blk] & bit_mask) != 0;
	}
	if ((field.ptr[beg_blk] & (~(detail::bit_pattern(range.beg) - 1))) != 0) return true;
	if ((field.ptr[end_blk] & (detail::bit_pattern(range.end) - 1)) != 0) return true;

	const int64 size = end_blk - beg_blk - 1;
	if (size > 0) {
		const uint8* p = (uint8*)(field.ptr + beg_blk + 1);
		const auto s = size * sizeof(Bitfield::ElementType);
		return !(p[0] == 0 && !memcmp(p, p + 1, s - 1));
	}

	return false;
}

inline bool all_bits_set(const Bitfield field) {
	const uint8* p = (uint8*)(field.ptr);
	const auto s = field.size_in_bytes();
	return (p[0] == 0xFF && !memcmp(p, p + 1, s - 1));
}

template <typename Int>
inline bool all_bits_set_in_range(const Bitfield field, Range<Int> range) {
	const auto beg_blk = detail::block(range.beg);
	const auto end_blk = detail::block(range.end);

	if (beg_blk == end_blk) {
		// All bits reside within the same Block
		const auto bit_mask = (detail::bit_pattern(range.beg) - 1) ^ (detail::bit_pattern(range.end) - 1);
		return (field.ptr[beg_blk] & bit_mask) == bit_mask;
	}

	if ((field.ptr[beg_blk] & (~(detail::bit_pattern(range.beg) - 1))) != (~(detail::bit_pattern(range.beg) - 1))) return false;
	if ((field.ptr[end_blk] & (detail::bit_pattern(range.end) - 1)) != detail::bit_pattern(range.end) - 1) return false;

	const int64 size = end_blk - beg_blk - 1;
	if (size > 0) {
		const uint8* p = (uint8*)(field.ptr + beg_blk + 1);
		const auto s = size * sizeof(Bitfield::ElementType);
		return (p[0] == 0xFF && !memcmp(p, p + 1, s - 1));
	}

	return true;
}

inline bool get_bit(const Bitfield field, int64 idx) {
    return (field.ptr[detail::block(idx)] & detail::bit_pattern(idx)) != 0U;
}

inline void set_bit(Bitfield field, int64 idx) {
    field.ptr[detail::block(idx)] |= detail::bit_pattern(idx);
}

inline void clear_bit(Bitfield field, int64 idx) {
	field.ptr[detail::block(idx)] &= ~detail::bit_pattern(idx);
}

inline bool invert_bit(Bitfield field, int64 idx) {
    return field.ptr[detail::block(idx)] ^= detail::bit_pattern(idx);
}

inline void and_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
	// @TODO: Vectorize
    for (int64 i = 0; i < detail::num_blocks(dst); i++) {
        dst.ptr[i] = src_a.ptr[i] & src_b.ptr[i];
    }
}

inline void or_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
	// @TODO: Vectorize
    for (int64 i = 0; i < detail::num_blocks(dst); i++) {
        dst.ptr[i] = src_a.ptr[i] | src_b.ptr[i];
    }
}

inline void xor_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
	ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
	// @TODO: Vectorize
	for (int64 i = 0; i < detail::num_blocks(dst); i++) {
		dst.ptr[i] = src_a.ptr[i] ^ src_b.ptr[i];
	}
}

template <typename T>
int64 extract_data_from_mask(T* RESTRICT out_data, const T* RESTRICT in_data, Bitfield mask) {
	constexpr uint32 bits_per_block = 32;
	STATIC_ASSERT(bits_per_block == (sizeof(Bitfield::ElementType) * 8), "Bits per block does not match ElementType");
	int64 out_count = 0;
	for (int64 i = 0; i < mask.size(); i += bits_per_block) {
		const auto block = mask.ptr[i / bits_per_block];
		if (block != 0U) {
			for (; i < mask.size(); i++) {
				if (bitfield::get_bit(mask, i)) {
					out_data[out_count] = in_data[i];
					out_count++;
				}
			}
		}
	}
	return out_count;
}

void print(const Bitfield field);

}  // namespace Bitfield
