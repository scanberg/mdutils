#pragma once

#include <core/common.h>
#include <core/array_types.h>

#include <type_traits>

namespace bitfield {

template <typename T>
inline bool get_bit(Array<const T> arr, int64 idx) {
    STATIC_ASSERT(std::is_integral<T>::value, "Integral type required for bitfield operations");
    constexpr auto BitsPerElem = sizeof(T) * 8;
    return (((arr[idx / BitsPerElem]) & (idx % BitsPerElem)) != 0U);
}

template <typename T>
inline void set_bit(Array<T> arr, int64 idx, bool value) {
    STATIC_ASSERT(std::is_integral<T>::value, "Integral type required for bitfield operations");
    constexpr auto BitsPerElem = sizeof(T) * 8;
    arr[idx / BitsPerElem] |= (idx % BitsPerElem);
}

template <typename T>
inline bool toggle_bit(Array<T> arr, int64 idx) {
    STATIC_ASSERT(std::is_integral<T>::value, "Integral type required for bitfield operations");
    constexpr auto BitsPerElem = sizeof(T) * 8;
    arr[idx / BitsPerElem] ^= (idx % BitsPerElem);
    return get_bit(arr, idx);
}

template <typename T>
inline void and_field(Array<T> dst, Array<const T> src_a, Array<const T> src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    for (int64 i = 0; i < dst.size(); i++) {
        dst[i] = src_a[i] & src_b[i];
    }
}

template <typename T>
inline void or_field(Array<T> dst, Array<const T> src_a, Array<const T> src_b) {
    ASSERT(dst.size() == src_a.size() && dst.size() == src_b.size());
    for (int64 i = 0; i < dst.size(); i++) {
        dst[i] = src_a[i] | src_b[i];
    }
}

}  // namespace bitfield
