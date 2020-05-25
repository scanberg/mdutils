#pragma once

#include <core/common.h>
#include <core/types.h>
#include <string.h>  // memset, memcpy
//#include <core/array_types.h>

struct Bitfield {
    using BlockType = u64;

    BlockType* block_ptr = nullptr;
    i64 bit_count = 0;

    constexpr i64 size() const { return bit_count; }
    constexpr i64 size_in_bytes() const { return ((bit_count + (sizeof(BlockType) * 8) - 1) / (sizeof(BlockType) * 8)) * sizeof(BlockType); }  // @NOTE: round up integer div.

    constexpr BlockType* data() { return block_ptr; }
    constexpr BlockType* begin() { return block_ptr; }
    constexpr BlockType* beg() { return block_ptr; }
    constexpr BlockType* end() { return block_ptr + bit_count; }

    constexpr const BlockType* data() const { return block_ptr; }
    constexpr const BlockType* begin() const { return block_ptr; }
    constexpr const BlockType* beg() const { return block_ptr; }
    constexpr const BlockType* end() const { return block_ptr + bit_count; }

    constexpr bool empty() const { return bit_count > 0; }
    constexpr bool operator[](const i64 bit_idx) const;

    constexpr operator bool() const { return block_ptr != nullptr; }
};

namespace bitfield {

namespace detail {
constexpr auto bits_per_block = sizeof(Bitfield::BlockType) * 8;

constexpr Bitfield::BlockType block_idx(i64 bit_idx) { return (Bitfield::BlockType)(bit_idx / bits_per_block); }
constexpr Bitfield::BlockType bit_pattern(i64 bit_idx) { return (Bitfield::BlockType)1 << (bit_idx % bits_per_block); }

// from https://yesteapea.wordpress.com/2013/03/03/counting-the-number-of-set-bits-in-an-integer/
constexpr i64 number_of_set_bits(u64 i) {
    i = i - ((i >> 1) & 0x5555555555555555);
    i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
    i = ((i + (i >> 4)) & 0x0F0F0F0F0F0F0F0F);
    return (i * (0x0101010101010101)) >> 56;
}

constexpr i64 num_blocks(Bitfield field) {
    return (field.size_in_bytes() + sizeof(Bitfield::BlockType) - 1) / sizeof(Bitfield::BlockType);  // @NOTE: round up integer div.
}

// @TODO: Use intrinsics here perhaps, this should be 1 instruction on modern cpus
constexpr i64 bit_scan_forward(Bitfield::BlockType mask) {
    if (mask == 0) return -1;
    i64 idx = 0;
    while ((mask & 1) == 0) {
        mask = mask >> 1;
        ++idx;
    }
    return idx;
}

constexpr i64 bit_scan_reverse(Bitfield::BlockType mask) {
    if (mask == 0) return -1;
    i64 idx = bits_per_block - 1;
    constexpr Bitfield::BlockType last_bit = ((Bitfield::BlockType)1 << (bits_per_block - 1));
    while ((mask & last_bit) == 0) {
        mask = mask << 1;
        --idx;
    }
    return idx;
}

}  // namespace detail

inline void free(Bitfield* field) {
    ASSERT(field);
    if (field->block_ptr) {
        ALIGNED_FREE(field->block_ptr);
        field->block_ptr = nullptr;
        field->bit_count = 0;
    }
}

inline void init(Bitfield* field, i64 num_bits) {
    ASSERT(field);
    free(field);
    field->bit_count = num_bits;
    const i64 alloc_size = detail::num_blocks(*field) * sizeof(Bitfield::BlockType);
    field->block_ptr = (Bitfield::BlockType*)ALIGNED_MALLOC(alloc_size, 16);  // @NOTE: Align to 16 byte to allow for aligned simd load/store
    memset(field->block_ptr, 0, alloc_size);
}

inline void init(Bitfield* field, Bitfield src) {
    ASSERT(field);
    free(field);
    if (!src) return;
    field->block_ptr = (Bitfield::BlockType*)ALIGNED_MALLOC(src.size_in_bytes(), 16);
    field->bit_count = src.bit_count;
    memcpy(field->block_ptr, src.block_ptr, src.size_in_bytes());
}

inline void copy(Bitfield dst, const Bitfield src) {
    ASSERT(dst.size() == src.size() && "Bitfield size did not match");
    memcpy(dst.block_ptr, src.block_ptr, dst.size_in_bytes());
}

inline void set_all(Bitfield field) { memset(field.block_ptr, 0xFF, field.size_in_bytes()); }

inline void clear_all(Bitfield field) { memset(field.block_ptr, 0, field.size_in_bytes()); }

constexpr void invert_all(Bitfield field) {
    // @TODO: Vectorize
    for (i64 i = 0; i < detail::num_blocks(field); i++) {
        field.block_ptr[i] = ~field.block_ptr[i];
    }
}

constexpr i64 number_of_bits_set(const Bitfield field) {
    const u32* ptr = (u32*)field.data();
    const i64 stride = sizeof(u32) * 8;
    const i64 size = field.size() / stride;
    const i64 rest = field.size() % stride;
    i64 count = 0;
    for (i64 i = 0; i < size; i++) {
        count += detail::number_of_set_bits(ptr[i]);
    }
    if (rest != 0) {
        const auto bits = ptr[size] & (detail::bit_pattern(rest) - 1);
        count += detail::number_of_set_bits(bits);
    }

    return count;
}

template <typename Int>
constexpr void set_range(Bitfield field, Range<Int> range) {
    const auto beg_blk = detail::block_idx(range.beg);
    const auto end_blk = detail::block_idx(range.end);
    const auto beg_mask = ~(detail::bit_pattern(range.beg) - 1);
    const auto end_mask = (detail::bit_pattern(range.end) - 1);

    if (beg_blk == end_blk) {
        // All bits reside within the same Block
        const auto mask = beg_mask & end_mask;
        field.block_ptr[beg_blk] |= mask;
        return;
    }

    field.block_ptr[beg_blk] |= beg_mask;
    // Set any bits within the inner range of blocks: beg_blk, [inner range], end_blk
    const i64 size = end_blk - beg_blk - 1;
    if (size > 0) {
        memset(field.block_ptr + beg_blk + 1, 0xFF, size * sizeof(Bitfield::BlockType));
    }
    field.block_ptr[end_blk] |= end_mask;
}

template <typename Int>
constexpr void clear_range(Bitfield field, Range<Int> range) {
    const auto beg_blk = detail::block_idx(range.beg);
    const auto end_blk = detail::block_idx(range.end);
    const auto beg_mask = ~(detail::bit_pattern(range.beg) - 1);
    const auto end_mask = (detail::bit_pattern(range.end) - 1);

    if (beg_blk == end_blk) {
        // All bits reside within the same Block
        const auto mask = beg_mask & end_mask;
        field.block_ptr[beg_blk] &= ~mask;
        return;
    }


    field.block_ptr[beg_blk] &= ~beg_mask;
    // Set any bits within the inner range of blocks: beg_blk, [inner range], end_blk
    const i64 size = end_blk - beg_blk - 1;
    if (size > 0) {
        memset(field.block_ptr + beg_blk + 1, 0x00, size * sizeof(Bitfield::BlockType));
    }
    field.block_ptr[end_blk] &= ~end_mask;
}

template <typename Int>
constexpr bool any_bit_set_in_range(const Bitfield field, Range<Int> range) {
    const auto beg_blk = detail::block_idx(range.beg);
    const auto end_blk = detail::block_idx(range.end);
    const auto beg_mask = ~(detail::bit_pattern(range.beg) - 1);
    const auto end_mask = (detail::bit_pattern(range.end) - 1);

    if (beg_blk == end_blk) {
        // All bits reside within the same Block
        const auto mask = beg_mask & end_mask;
        return (field.block_ptr[beg_blk] & mask) != 0;
    }

    // Mask out and explicitly check beg and end blocks
    if ((field.block_ptr[beg_blk] & beg_mask) != 0) return true;
    if ((field.block_ptr[end_blk] & end_mask) != 0) return true;

    // memcmp rest
    const i64 size = end_blk - beg_blk - 1;
    if (size > 0) {
        const u8* p = (u8*)(field.block_ptr + beg_blk + 1);
        const auto s = size * sizeof(Bitfield::BlockType);
        return !(p[0] == 0 && !memcmp(p, p + 1, s - 1));
    }

    return false;
}

constexpr bool any_bit_set(const Bitfield field) {
    constexpr Bitfield::BlockType beg_blk_idx = 0;
    const     Bitfield::BlockType end_blk_idx = detail::block_idx(field.bit_count);

    if (beg_blk_idx == end_blk_idx) {
        // All bits reside within the same Block
        const auto mask = detail::bit_pattern(field.bit_count) - 1;
        return (field.block_ptr[beg_blk_idx] & mask) != 0;
    }

    for (u64 i = 0; i < end_blk_idx - 1; i++) {
        if (field.block_ptr[i] != 0) return true;
    }
    if ((field.block_ptr[end_blk_idx] & (detail::bit_pattern(field.bit_count) - 1)) != 0) return true;
    return false;
}

constexpr bool all_bits_set(const Bitfield field) {
    const auto beg_blk = 0;
    const auto end_blk = detail::block_idx(field.size());
    const auto end_mask = (detail::bit_pattern(field.size()) - 1);

    if (beg_blk == end_blk) {
        return (field.block_ptr[beg_blk] & end_mask) == end_mask;
    }

    for (u64 i = 0; i < end_blk - 1; i++) {
        if (~field.block_ptr[i] != 0) return false;
    }
    return (field.block_ptr[end_blk] & end_mask) == end_mask;
}

template <typename Int>
constexpr bool all_bits_set_in_range(const Bitfield field, Range<Int> range) {
    const auto beg_blk = detail::block_idx(range.beg);
    const auto end_blk = detail::block_idx(range.end);
    const auto beg_mask = ~(detail::bit_pattern(range.beg) - 1);
    const auto end_mask = (detail::bit_pattern(range.end) - 1);

    if (beg_blk == end_blk) {
        // All bits reside within the same Block
        const auto mask = beg_mask & end_mask;
        return (field.block_ptr[beg_blk] & mask) == mask;
    }

    // Mask out and explicitly check beg and end blocks
    if ((field.block_ptr[beg_blk] & beg_mask) != beg_mask) return false;
    if ((field.block_ptr[end_blk] & end_mask) != end_mask) return false;

    const i64 size = end_blk - beg_blk - 1;
    if (size > 0) {
        const u8* p = (u8*)(field.block_ptr + beg_blk + 1);
        const auto s = size * sizeof(Bitfield::BlockType);
        return (p[0] == 0xFF && !memcmp(p, p + 1, s - 1));
    }

    return true;
}

// Finds the next bit set in the field, beggining with a supplied offset which is included in search
constexpr i64 find_next_bit_set(const Bitfield field, i64 offset = 0) {
    if (offset >= field.size()) return -1;

    Bitfield::BlockType blk_idx = (offset / detail::bits_per_block);
    const Bitfield::BlockType num_blocks = detail::num_blocks(field);

    // Check first block explicitly with proper mask
    auto mask = field.block_ptr[blk_idx] & (~(detail::bit_pattern(offset) - 1));
    if (mask != 0) {
        return detail::bit_scan_forward(mask);
    }

    offset += detail::bits_per_block;
    for (++blk_idx; blk_idx < num_blocks - 1; blk_idx++) {
        mask = field.block_ptr[blk_idx];
        if (mask != 0) {
            return offset + detail::bit_scan_forward(mask);
        }
        offset += detail::bits_per_block;
    }

    mask = field.block_ptr[blk_idx] & (detail::bit_pattern(field.bit_count) - 1);
    if (mask != 0) {
        return offset + detail::bit_scan_forward(mask);
    }

    return -1;
}

constexpr i64 find_first_bit_set(const Bitfield field) {
    for (i64 blk_idx = 0; blk_idx < detail::num_blocks(field) - 1; blk_idx++) {
        const auto blk = field.block_ptr[blk_idx];
        if (blk != 0) {
            return blk_idx * detail::bits_per_block + detail::bit_scan_forward(blk);
        }
    }

    const auto last_blk_idx = detail::num_blocks(field) - 1;
    const auto last_blk = field.block_ptr[last_blk_idx] & (detail::bit_pattern(field.bit_count) - 1);
    if (last_blk != 0) {
        return last_blk_idx * detail::bits_per_block + detail::bit_scan_forward(last_blk);
    }

    return -1;
}

constexpr i64 find_last_bit_set(const Bitfield field) {
    const auto last_blk_idx = detail::num_blocks(field) - 1;
    const auto last_blk = field.block_ptr[last_blk_idx] & (detail::bit_pattern(field.bit_count) - 1);
    if (last_blk != 0) {
        return last_blk_idx * detail::bits_per_block + detail::bit_scan_reverse(last_blk);
    }

    for (i64 blk_idx = detail::num_blocks(field) - 2; blk_idx >= 0; blk_idx--) {
        const auto blk = field.block_ptr[blk_idx];
        if (blk != 0) {
            return blk_idx * detail::bits_per_block + detail::bit_scan_reverse(blk);
        }
    }

    return -1;
}

constexpr bool get_bit(const Bitfield field, i64 idx) { return (field.block_ptr[detail::block_idx(idx)] & detail::bit_pattern(idx)) != 0U; }

constexpr void set_bit(Bitfield field, i64 idx) { field.block_ptr[detail::block_idx(idx)] |= detail::bit_pattern(idx); }

constexpr void clear_bit(Bitfield field, i64 idx) { field.block_ptr[detail::block_idx(idx)] &= ~detail::bit_pattern(idx); }

constexpr bool invert_bit(Bitfield field, i64 idx) { return field.block_ptr[detail::block_idx(idx)] ^= detail::bit_pattern(idx); }

inline void and_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    // @TODO: Vectorize
    for (i64 i = 0; i < detail::num_blocks(dst); i++) {
        dst.block_ptr[i] = src_a.block_ptr[i] & src_b.block_ptr[i];
    }
}

inline void and_not_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    // @TODO: Vectorize
    for (i64 i = 0; i < detail::num_blocks(dst); i++) {
        dst.block_ptr[i] = src_a.block_ptr[i] & ~src_b.block_ptr[i];
    }
}

inline void or_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    // @TODO: Vectorize
    for (i64 i = 0; i < detail::num_blocks(dst); i++) {
        dst.block_ptr[i] = src_a.block_ptr[i] | src_b.block_ptr[i];
    }
}

inline void or_not_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    // @TODO: Vectorize
    for (i64 i = 0; i < detail::num_blocks(dst); i++) {
        dst.block_ptr[i] = src_a.block_ptr[i] | ~src_b.block_ptr[i];
    }
}

inline void xor_field(Bitfield dst, const Bitfield src_a, const Bitfield src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    // @TODO: Vectorize
    for (i64 i = 0; i < detail::num_blocks(dst); i++) {
        dst.block_ptr[i] = src_a.block_ptr[i] ^ src_b.block_ptr[i];
    }
}

template <typename T>
int64_t gather_masked(T* dst_data, const T* src_data, Bitfield mask, int64_t offset = 0) {
    int64_t dst_count = 0;
    const int64_t last_blk = detail::num_blocks(mask) - 1;
    for (int64_t blk_idx = 0; blk_idx < last_blk; ++blk_idx) {
        const auto blk = mask.block_ptr[blk_idx];
        if (blk == 0) continue;

        for (uint64_t i = 0; i < detail::bits_per_block; ++i) {
            if (blk & ((Bitfield::BlockType)1 << i)) {
                dst_data[dst_count] = src_data[offset + blk_idx * detail::bits_per_block + i];
                ++dst_count;
            }
        }
    }

    // @NOTE: Remainder
    const uint64_t blk_idx = last_blk;
    const uint64_t remainder = (uint64_t)mask.size() - blk_idx * detail::bits_per_block;
    const uint64_t blk = mask.block_ptr[blk_idx];
    if (blk) {
        for (uint64_t i = 0; i < remainder; ++i) {
            const uint64_t bit = (Bitfield::BlockType)1 << i;
            if (blk & bit) {
                dst_data[dst_count] = src_data[offset + blk_idx * detail::bits_per_block + i];
                ++dst_count;
            }
        }
    }

    return dst_count;
}

template <typename T>
int64_t scatter_masked(T* dst_data, const T* src_data, Bitfield mask, int64_t offset = 0) {
    int64_t src_count = 0;
    const int64_t last_blk = detail::num_blocks(mask) - 1;
    for (int64_t blk_idx = 0; blk_idx < last_blk; ++blk_idx) {
        const auto blk = mask.block_ptr[blk_idx];
        if (blk == 0) continue;

        for (uint64_t i = 0; i < detail::bits_per_block; ++i) {
            if (blk & ((Bitfield::BlockType)1 << i)) {
                dst_data[offset + blk_idx * detail::bits_per_block + i] = src_data[src_count];
                ++src_count;
            }
        }
    }

    // @NOTE: Remainder
    const uint64_t blk_idx = last_blk;
    const uint64_t remainder = (uint64_t)mask.size() - blk_idx * detail::bits_per_block;
    const uint64_t blk = mask.block_ptr[blk_idx];
    if (blk) {
        for (uint64_t i = 0; i < remainder; ++i) {
            const uint64_t bit = (Bitfield::BlockType)1 << i;
            if (blk & bit) {
                dst_data[offset + blk_idx * detail::bits_per_block + i] = src_data[src_count];
                ++src_count;
            }
        }
    }

    return src_count;
}

template <typename Func>
void for_each_bit_set(Bitfield mask, Func f, int64_t offset = 0) {
    ASSERT(mask);
    const int64_t last_blk = detail::num_blocks(mask) - 1;
    for (int64_t blk_idx = 0; blk_idx < last_blk; ++blk_idx) {
        const auto blk = mask.block_ptr[blk_idx];
        if (blk == 0) continue;

        for (uint64_t i = 0; i < detail::bits_per_block; ++i) {
            const uint64_t bit = (Bitfield::BlockType)1 << i;
            if (blk & bit) {
                const int64_t idx = offset + blk_idx * detail::bits_per_block + i;
                f(idx);
            }
        }
    }

    // @NOTE: Remainder
    const uint64_t blk_idx = last_blk;
    const uint64_t remainder = (uint64_t)mask.size() - blk_idx * detail::bits_per_block;
    const uint64_t blk = mask.block_ptr[blk_idx];
    if (blk) {
        for (uint64_t i = 0; i < remainder; ++i) {
            const uint64_t bit = (Bitfield::BlockType)1 << i;
            if (blk & bit) {
                const int64_t idx = offset + blk_idx * detail::bits_per_block + i;
                f(idx);
            }
        }
    }
}

void print(const Bitfield field);

}  // namespace bitfield

constexpr bool Bitfield::operator[](const i64 bit_idx) const { return block_ptr[bit_idx / bitfield::detail::bits_per_block] & ((Bitfield::BlockType)1 << (bit_idx % bitfield::detail::bits_per_block)); }
