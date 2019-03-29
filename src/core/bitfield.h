#pragma once

#include <core/common.h>
#include <core/types.h>
#include <core/array_types.h>

struct Bitfield {
	typedef uint64 ElementType;

	ElementType* ptr = nullptr;
	int64 count = 0;

	int64 size() const { return count; }
	int64 size_in_bytes() const { return count * sizeof(ElementType); }
	
	ElementType* data() { return ptr; }
	ElementType* beg() { return ptr; }
	ElementType* end() { return ptr + count; }

	operator bool() const {
		return ptr != nullptr && count != 0;
	}
};
//typedef Array<uint64> Bitfield;

namespace bitfield {

namespace detail {
constexpr auto block_bits = sizeof(Bitfield::ElementType) * 8;

inline Bitfield::ElementType block(int64 idx) {
	return idx / block_bits;
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
	return field.size() / block_bits;
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

inline void init(Bitfield* field, int64 num_bits, bool value = false) {
    ASSERT(field);
    free(field);
    const auto num_blocks = num_bits / detail::block_bits + 1;
    field->ptr = (uint64*)ALIGNED_MALLOC(num_blocks * sizeof(Bitfield::ElementType) + 16, 16);   // @NOTE: We pad and align to 16 byte to allow for aligned simd load/store
	field->count = num_bits;
	memset(field->ptr, (int)value, field->size_in_bytes());
}

inline void init(Bitfield* field, Bitfield src) {
	ASSERT(field);
	free(field);
	if (!src) return;
	field->ptr = (uint64*)ALIGNED_MALLOC(src.size_in_bytes(), 16);
	field->count = src.count;
	memcpy(field->ptr, src.ptr, src.size_in_bytes());
}

inline void set_all(Bitfield field) {
    memset(field.ptr, 1, field.size_in_bytes());
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

	// Set bits inside beg_blk
	field.ptr[beg_blk] |= ~(detail::bit_pattern(range.beg) - 1);

	// Set bits inside end_blk
	field.ptr[beg_blk] |= (detail::bit_pattern(range.end) - 1);

	// Set any bits within the inner range of blocks: beg_blk, [inner range], end_blk
	const int64 size = end_blk - beg_blk - 2;
	if (size > 0) {
		memset(field.ptr + beg_blk, 1, size * sizeof(Bitfield::ElementType));
	}
}

inline bool any_bit_set(Bitfield field) {
	const auto size = field.size_in_bytes();
	const uint8* buf = (uint8*)field.data();
	return !(buf[0] == 0 && !memcmp(buf, buf + 1, size - 1));
}

template <typename Int>
inline bool any_bit_set_in_range(Bitfield field, Range<Int> range) {
	const auto beg_blk = detail::block(range.beg);
	const auto end_blk = detail::block(range.end);

	if (beg_blk == end_blk) {
		// All bits reside within the same Block
		const auto bits = (detail::bit_pattern(range.beg) - 1) ^ (bit_pattern(range.end) - 1);
		return field.ptr[beg_blk] & bits != 0;
	}
	if (field.ptr[beg_blk] & (~(detail::bit_pattern(range.beg) - 1)) != 0) return true;
	if (field.ptr[end_blk] & (detail::bit_pattern(range.end) - 1) != 0) return true;

	const int64 size = end_blk - beg_blk - 2;
	if (size > 0) {
		const uint8* p = (uint8*)field.ptr[beg_blk + 1];
		const auto s = size * sizeof(Bitfield::ElementType);
		return !(p[0] == 0 && !memcmp(p, p + 1, s - 1));
	}

	return false;
}

inline int64 number_of_bits_set(Bitfield field) {
	const uint32* ptr = (uint32*)field.data();
	const int64 num_blocks_32 = field.size() / 32 + 1;
	int64 count = 0;
	for (int64 i = 0; i < num_blocks_32; i++) {
		count += detail::number_of_set_bits(ptr[i]);
	}
	return count;
}

inline bool get_bit(Bitfield field, int64 idx) {
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

void print(Bitfield field);

}  // namespace Bitfield
