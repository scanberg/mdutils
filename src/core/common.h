#pragma once

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>

#ifndef NDEBUG
#define DEBUG
/*
inline bool _assert(const char* file, const char* func, int line, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "ASSERTION FAILED!\n%s:%s:%i\n", file, func, line);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    assert(false);
    return false;
}
inline bool _assert(const char* file, const char* func, int line) { _assert(file, func, line, ""); }
*/
#define ASSERT(COND) assert(COND)
#else
#define ASSERT(...) \
    {}
#endif
#define STATIC_ASSERT(cond, msg) static_assert(cond, msg)

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

#define KILOBYTES(x) (x << 10)
#define MEGABYTES(x) (KILOBYTES(x) << 10)
#define GIGABYTES(x) (MEGABYTES(x) << 10)

#define BIT(bit) (1U << bit)
#define SET_BIT(field, bit) (field |= (1U << bit))
#define CLEAR_BIT(field, bit) (field &= ~(1U << bit))
#define CHECK_BIT(field, bit) ((field & (1U << bit)) != 0U)

#define UNUSED(x) (void)(x)

#ifndef MALLOC
#define MALLOC(size) malloc(size)
#if _MSC_VER && !__INTEL_COMPILER
#include <malloc.h>
#define ALIGNED_MALLOC(size, alignment) _aligned_malloc(size, alignment)
#define ALIGNED_FREE(addr) _aligned_free(addr)
#else
#ifdef __GNUC__
#include <mm_malloc.h>
#endif
#define ALIGNED_MALLOC(size, alignment) _mm_malloc(size, alignment)
#define ALIGNED_FREE(addr) _mm_free(addr)
#endif
#define REALLOC(ptr, new_size) realloc(ptr, new_size)
#define CALLOC(num_items, size_of_item) calloc(num_items, size_of_item)
#define FREE(addr) free(addr)
#endif

// This is bogus atm and should be implemented properly
#ifndef TMP_MALLOC
#define TMP_MALLOC(size) malloc(size)
#if _MSC_VER && !__INTEL_COMPILER
#include <malloc.h>
#define TMP_ALIGNED_MALLOC(size, alignment) _aligned_malloc(size, alignment)
#define TMP_ALIGNED_FREE(addr) _aligned_free(addr)
#else
#define TMP_ALIGNED_MALLOC(size, alignment) _mm_malloc(size, alignment)
#define TMP_ALIGNED_FREE(addr) _mm_free(addr)
#endif
#define TMP_REALLOC(ptr, new_size) realloc(ptr, new_size)
#define TMP_CALLOC(num_items, size_of_item) calloc(num_items, size_of_item)
#define TMP_FREE(addr) free(addr)
#endif

#ifdef NEW
#undef NEW
#endif

#ifdef PLACEMENT_NEW
#undef PLACEMENT_NEW
#endif

#ifdef DELETE
#undef DELETE
#endif

#ifdef TMP_NEW
#undef TMP_NEW
#endif

#ifdef TMP_DELETE
#undef TMP_DELETE
#endif

// Blatantly stolen from ImGui (thanks Omar!)
struct NewDummy {};
inline void* operator new(size_t, NewDummy, void* ptr) { return ptr; }
inline void  operator delete(void*, NewDummy, void*)   {} // This is only required so we can use the symetrical new()
#define PLACEMENT_NEW(_PTR)					new(NewDummy(), _PTR)
#define NEW(_TYPE)							new(NewDummy(), MALLOC(sizeof(_TYPE))) _TYPE
#define TMP_NEW(_TYPE)						new(NewDummy(), TMP_MALLOC(sizeof(_TYPE))) _TYPE
template<typename T> void DELETE(T* p)		{ if (p) { p->~T(); FREE(p); } }
template<typename T> void TMP_DELETE(T* p)	{ if (p) { p->~T(); TMP_FREE(p); } }

#define RESTRICT __restrict

inline void* get_next_aligned_adress(void* mem, size_t align) {
    const size_t addr = (size_t)mem;
    return (void*)((addr + (align - 1)) & (~align + 1));
}

#define IS_ALIGNED(ptr, alignment) (((uintptr_t)ptr % alignment) == 0)

// Implementation of 'defer' in c++.
// From here https://pastebin.com/suTkpYp4
// With some slight modifications

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
