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

inline uint32_t clz(uint32_t v) {
    return __lzcnt(v);
}

inline uint64_t clz(uint64_t v) {
    return __lzcnt64(v);
}

#endif

inline int find_first_zero_byte(uint32_t x) {
    uint32_t y = (x & 0x7F7F7F7F) + 0x7F7F7F7F;
    y = ~(y | x | 0x7F7F7F7F);
    return (31U - clz(y)) >> 3;
}

inline int find_first_zero_byte(uint64_t x) {
    //uint64_t y = (x & 0x7F7F7F7F7F7F7F7F) + 0x7F7F7F7F7F7F7F7F;
    //y = ~(y | x | 0x7F7F7F7F7F7F7F7F);
    //return (63LLU - clz(y)) >> 3LLU;
    uint64_t y = (x - 0x0101010101010101) & ~x & 0x8080808080808080;
    return clz(y) >> 3;
}