#pragma once

#include <core/common.h>
#include <core/array_types.h>

typedef Array<uint64> Bitfield;

namespace bitfield {

constexpr auto block_bits = sizeof(Bitfield::ElementType) * 8;

inline Bitfield::ElementType block(int64 idx) {
	return idx / block_bits;
}

inline Bitfield::ElementType bit(int64 idx) {
	return (Bitfield::ElementType)1 << (idx % block_bits);
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
    const auto num_blocks = num_bits / block_bits + 1;
    field->ptr = (uint64*)ALIGNED_MALLOC(num_blocks * sizeof(Bitfield::ElementType) + 16, 16);   // @NOTE: We pad and align to 16 byte to aid simd operations
	field->count = num_bits / block_bits;
	memset(field->ptr, (int)value, field->size_in_bytes());
}

inline void set_all(Bitfield field, bool value) {
    memset(field.ptr, (int)value, field.size_in_bytes());
}

inline void clear(Bitfield field) {
    memset(field.ptr, 0, field.size_in_bytes());
}

inline void invert(Bitfield field) {
    // @TODO: Vectorize 
    for (int64 i = 0; i < field.size(); i++) {
        field[i] = ~field[i];
    }
}

inline bool get_bit(Bitfield field, int64 idx) {
    return (field[block(idx)] & bit(idx)) != 0U;
}

inline void set_bit(Bitfield field, int64 idx) {
    field[block(idx)] |= bit(idx);
}

inline void clear_bit(Bitfield field, int64 idx) {
	field[block(idx)] &= ~bit(idx);
}

inline bool toggle_bit(Bitfield field, int64 idx) {
    return field[block(idx)] ^= bit(idx);
}

inline void and_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
	// @TODO: Vectorize
    for (int64 i = 0; i < dst.size(); i++) {
        dst[i] = src_a[i] & src_b[i];
    }
}

inline void or_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
	// @TODO: Vectorize
    for (int64 i = 0; i < dst.size(); i++) {
        dst[i] = src_a[i] | src_b[i];
    }
}

void print(Bitfield field);

}  // namespace Bitfield
