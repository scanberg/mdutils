#pragma once

#include <core/common.h>
#include <core/array_types.h>

#include <type_traits>

typedef Array<uint64> BitField;

namespace bitfield {

constexpr auto block_bits = sizeof(BitField::ElementType);
constexpr auto simd_block_bits = 128;

inline void free_bitfield(BitField* field) {
    ASSERT(field);
    if (field->ptr) {
        ALIGNED_FREE(field->ptr);
        field->ptr = nullptr;
        field->count = 0;
    }
}

inline void init_bitfield(BitField* field, int64 num_bits) {
    ASSERT(field);
    free_bitfield(field);
    constexpr auto num_blocks = num_bits / simd_block_bits + 1;
    field->ptr = ALIGNED_MALLOC(num_blocks * simd_block_bytes, simd_block_bytes);   // @NOTE: We pad and align to 16 byte to aid simd operations
}

inline void set_all(BitField field, bool value) {
    memset(field.ptr, (int)value, field.size_in_bytes());
}

inline void clear(BitField field) {
    memset(field.ptr, 0, field.size_in_bytes());
}

inline void invert(BitField field) {
    // @TODO: Vectorize 
    for (int64 i = 0; i < field.size(); i++) {
        field[i] = ~field[i];
    }
}

inline bool get_bit(BitField arr, int64 idx) {
    return (((arr[idx / block_bits]) & (idx % block_bits)) != 0U);
}

inline void set_bit(BitField field, int64 idx, bool value) {
    arr[idx / block_bits] |= (idx % block_bits);
}

inline bool toggle_bit(BitField field, int64 idx) {
    constexpr auto BitsPerElem = sizeof(T) * 8;
    arr[idx / BitsPerElem] ^= (idx % BitsPerElem);
    return get_bit(arr, idx);
}

inline void and_field(BitField dst, const BitField src_a, const BitField src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    for (int64 i = 0; i < dst.size(); i++) {
        dst[i] = src_a[i] & src_b[i];
    }
}

inline void or_field(BitField dst, const BitField src_a, const BitField src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    for (int64 i = 0; i < dst.size(); i++) {
        dst[i] = src_a[i] | src_b[i];
    }
}

}  // namespace bitfield
