#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>

#include <stdlib.h>

#ifndef NDEBUG
#define DEBUG
inline void _assert(const char* file, const char* func, int line, bool cond, const char* fmt, ...) {
    if (!cond) {
        va_list ap;
        va_start(ap, fmt);
        fprintf(stderr, "ASSERTION FAILED!\n%s:%s:%i\n", file, func, line);
        vfprintf(stderr, fmt, ap);
        va_end(ap);
        assert(false);
    }
}
inline void _assert(const char* file, const char* func, int line, bool cond) { _assert(file, func, line, cond, ""); }

#define ASSERT(...) _assert(__FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)
#else
#define ASSERT(...) \
    {}
#endif
#define STATIC_ASSERT(cond, msg) static_assert(cond, msg)

#define KILOBYTES(x) (x << 10)
#define MEGABYTES(x) (KILOBYTES(x) << 10)
#define GIGABYTES(x) (MEGABYTES(x) << 10)

#define BIT(bit) (1U << bit)
#define SET_BIT(field, bit) (field |= (1U << bit))
#define CLEAR_BIT(field, bit) (field &= ~(1U << bit))
#define CHECK_BIT(field, bit) ((field & (1U << bit)) != 0U)

#define UNUSED(x) (void)(x)

#ifndef MALLOC
  #define MALLOC(x) malloc(x)
#if _MSC_VER && !__INTEL_COMPILER
#include <malloc.h>
  #define ALIGNED_MALLOC(x, y) _aligned_malloc(x, y)
  #define ALIGNED_FREE(x) _aligned_free(x)
#else
  #define ALIGNED_MALLOC(x, y) aligned_alloc(x, y)
  #define ALIGNED_FREE(x) free(x)
#endif
  #define REALLOC(x, y) realloc(x, y)
  #define CALLOC(x, y) calloc(x, y)
  #define FREE(x) free(x)
#endif

// This is bogus atm and should be implemented properly
#ifndef TMP_MALLOC
#define TMP_MALLOC(x) malloc(x)
#if _MSC_VER && !__INTEL_COMPILER
#include <malloc.h>
  #define TMP_ALIGNED_MALLOC(x, y) _aligned_malloc(x, y)
  #define TMP_ALIGNED_FREE(x) _aligned_free(x)
#else
  #define TMP_ALIGNED_MALLOC(x, y) aligned_alloc(x, y)
  #define ALIGNED_FREE(x) free(x)
#endif
#define TMP_REALLOC(x, y) realloc(x, y)
#define TMP_CALLOC(x, y) calloc(x, y)
#define TMP_FREE(x) free(x)
#endif

// implementation of 'defer' in c++.
// from here https://pastebin.com/suTkpYp4

#define CONCAT_INTERNAL(x, y) x##y
#define CONCAT(x, y) CONCAT_INTERNAL(x, y)

struct ExitScopeHelp {
    template <typename T>
    struct ExitScope {
        T lambda;
        ExitScope(T lambda) : lambda(lambda) {}
        ~ExitScope() { lambda(); }
        ExitScope& operator=(const ExitScope&) = delete;
    };

    template <typename T>
    ExitScope<T> operator+(T t) {
        return t;
    }
};

#define defer [[maybe_unused]] const auto& CONCAT(defer__, __LINE__) = ExitScopeHelp() + [&]()
