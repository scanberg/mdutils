#pragma once

#include "platform.h"

#if (COMPILER_CLANG || COMPILER_GCC)
inline uint32_t clz(uint32_t v) {
    return __builtin_clz(v);
}

inline uint64_t clz(uint64_t v) {
    return __builtin_clzll(v);
}
#endif

#if COMPILER_MSVC
#include <intrin.h>
#include <stdint.h>

// count bits
inline uint32_t popcnt(uint32_t v) {
    return __popcnt(v);
}

// count bits 64-bit
inline uint64_t popcnt(uint64_t v) {
    return __popcnt64(v);
}

inline uint32_t clz(uint32_t x) {
    return __lzcnt(x);
}

inline uint64_t clz(uint64_t x) {
    return __lzcnt64(x);
}

inline uint32_t ctz(uint32_t x) {
    return popcnt(~x & (x-1));
}

inline uint64_t ctz(uint64_t x) {
    return popcnt(~x & (x-1LLU));
}

#endif

inline uint32_t find_first_zero_byte(uint32_t x) {
    uint32_t y = (x - 0x01010101) & ~x & 0x80808080;
    return ctz(y) >> 3;
}

inline uint64_t find_first_zero_byte(uint64_t x) {
    uint64_t y = (x - 0x0101010101010101) & ~x & 0x8080808080808080;
    return ctz(y) >> 3;
}