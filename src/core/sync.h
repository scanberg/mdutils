#pragma once

#include <core/platform.h>
#include <stdint.h>

// Synchronization primitives and operations

// Atomic operations
inline int32_t atomic_fetch(volatile int32_t* _ptr);
inline int32_t atomic_fetch_and_set(volatile int32_t* _ptr, int32_t _val);
inline int32_t atomic_fetch_and_add(volatile int32_t* _ptr, int32_t _addend);
inline int32_t atomic_compare_and_swap(volatile int32_t* _ptr, int32_t _old, int32_t _new);
inline int32_t atomic_test_and_add(volatile int32_t* _ptr, int32_t _test, int32_t _addend);

inline uint32_t atomic_fetch(volatile uint32_t* _ptr);
inline uint32_t atomic_fetch_and_set(volatile uint32_t* _ptr, uint32_t _val);
inline uint32_t atomic_fetch_and_add(volatile uint32_t* _ptr, uint32_t _addend);
inline uint32_t atomic_compare_and_swap(volatile uint32_t* _ptr, uint32_t _old, uint32_t _new);
inline uint32_t atomic_test_and_add(volatile uint32_t* _ptr, uint32_t _test, uint32_t _addend);

inline int64_t atomic_fetch(volatile int64_t* _ptr);
inline int64_t atomic_fetch_and_set(volatile int64_t* _ptr, int64_t _val);
inline int64_t atomic_fetch_and_add(volatile int64_t* _ptr, int64_t _addend);
inline int64_t atomic_compare_and_swap(volatile int64_t* _ptr, int64_t _old, int64_t _new);
inline int64_t atomic_test_and_add(volatile int64_t* _ptr, int64_t _test, int64_t _addend);

inline uint64_t atomic_fetch(volatile uint64_t* _ptr);
inline uint64_t atomic_fetch_and_set(volatile uint64_t* _ptr, uint64_t _val);
inline uint64_t atomic_fetch_and_add(volatile uint64_t* _ptr, uint64_t _addend);
inline uint64_t atomic_compare_and_swap(volatile uint64_t* _ptr, uint64_t _old, uint64_t _new);
inline uint64_t atomic_test_and_add(volatile uint64_t* _ptr, uint64_t _test, uint64_t _addend);

inline void* atomic_exchange_ptr(void** _ptr, void* _new);

// TODO: Mutex, Semaphore, Conditional


// IMPLEMENTATION

#if PLATFORM_WINDOWS
#include <intrin.h>

inline int32_t atomic_fetch(volatile int32_t* _ptr) { return (int32_t)_InterlockedCompareExchange((volatile long*)_ptr, 0, 0); }
inline int32_t atomic_fetch_and_set(volatile int32_t* _ptr, int32_t _val) { return (int32_t)_InterlockedExchange((volatile long*)_ptr, _val); }
inline int32_t atomic_fetch_and_add(volatile int32_t* _ptr, int32_t _addend) { return (int32_t)_InterlockedExchangeAdd((volatile long*)_ptr, _addend); }
inline int32_t atomic_compare_and_swap(volatile int32_t* _ptr, int32_t _old, int32_t _new) {
    return (int32_t)_InterlockedCompareExchange((volatile long*)_ptr, _old, _new);
}

inline uint32_t atomic_fetch(volatile uint32_t* _ptr) { return (uint32_t)_InterlockedCompareExchange((volatile long*)_ptr, 0, 0); }
inline uint32_t atomic_fetch_and_set(volatile uint32_t* _ptr, uint32_t _val) { return (uint32_t)_InterlockedExchange((volatile long*)_ptr, _val); }
inline uint32_t atomic_fetch_and_add(volatile uint32_t* _ptr, uint32_t _addend) { return (uint32_t)_InterlockedExchangeAdd((volatile long*)_ptr, _addend); }
inline uint32_t atomic_compare_and_swap(volatile uint32_t* _ptr, uint32_t _old, uint32_t _new) {
    return (uint32_t)_InterlockedCompareExchange((volatile long*)_ptr, _old, _new);
}

inline int64_t atomic_fetch(volatile int64_t* _ptr) { return (int64_t)_InterlockedCompareExchange64((volatile long long*)_ptr, 0, 0); }
inline int64_t atomic_fetch_and_set(volatile int64_t* _ptr, int64_t _val) { return (int64_t)_InterlockedExchange64((volatile long long*)_ptr, _val); }
inline int64_t atomic_fetch_and_add(volatile int64_t* _ptr, int64_t _addend) { return (int64_t)_InterlockedExchangeAdd64((volatile long long*)_ptr, _addend); }
inline int64_t atomic_compare_and_swap(volatile int64_t* _ptr, int64_t _old, int64_t _new) {
    return (int64_t)_InterlockedCompareExchange64((volatile long long*)_ptr, _old, _new);
}

inline uint64_t atomic_fetch(volatile uint64_t* _ptr) { return (uint64_t)_InterlockedCompareExchange64((volatile long long*)_ptr, 0, 0); }
inline uint64_t atomic_fetch_and_set(volatile uint64_t* _ptr, uint64_t _val) { return (uint64_t)_InterlockedExchange64((volatile long long*)_ptr, _val); }
inline uint64_t atomic_fetch_and_add(volatile uint64_t* _ptr, uint64_t _addend) { return (uint64_t)_InterlockedExchangeAdd64((volatile long long*)_ptr, _addend); }
inline uint64_t atomic_compare_and_swap(volatile uint64_t* _ptr, uint64_t _old, uint64_t _new) {
    return (uint64_t)_InterlockedCompareExchange64((volatile long long*)_ptr, _old, _new);
}

inline void* atomic_exchange_ptr(void** _ptr, void* _new) { return _InterlockedExchangePointer(_ptr, _new); }

#elif PLATFORM_LINUX || PLATFORM_OSX
// ASSUME POSIX

inline int32_t atomic_fetch(volatile int32_t* _ptr) { return (int32_t)__sync_fetch_and_add(_ptr, 0); }
inline int32_t atomic_fetch_and_set(volatile int32_t* _ptr, int32_t _val) {
    int32_t result = (int32_t)__sync_lock_test_and_set(_ptr, _val);
    __sync_lock_release(_ptr);
    return result;
}
inline int32_t atomic_fetch_and_add(volatile int32_t* _ptr, int32_t _addend) { return (int32_t)__sync_fetch_and_add(_ptr, _addend); }
inline int32_t atomic_compare_and_swap(volatile int32_t* _ptr, int32_t _old, int32_t _new) {
    return (int32_t)___sync_val_compare_and_swap(_ptr, _old, _new);
}

inline uint32_t atomic_fetch(volatile uint32_t* _ptr) { return (uint32_t)__sync_fetch_and_add(_ptr, 0); }
inline uint32_t atomic_fetch_and_set(volatile uint32_t* _ptr, uint32_t _val) {
    uint32_t result = (uint32_t)__sync_lock_test_and_set(_ptr, _val);
    __sync_lock_release(_ptr);
    return result;
}
inline uint32_t atomic_fetch_and_add(volatile uint32_t* _ptr, uint32_t _addend) { return (uint32_t)__sync_fetch_and_add(_ptr, _addend); }
inline uint32_t atomic_compare_and_swap(volatile uint32_t* _ptr, uint32_t _old, uint32_t _new) {
    return (uint32_t)___sync_val_compare_and_swap(_ptr, _old, _new);
}

inline int64_t atomic_fetch(volatile int64_t* _ptr) { return (int64_t)__sync_fetch_and_add(_ptr, 0); }
inline int64_t atomic_fetch_and_set(volatile int64_t* _ptr, int64_t _val) {
    int64_t result = (int64_t)__sync_lock_test_and_set(_ptr, _val);
    __sync_lock_release(_ptr);
    return result;
}
inline int64_t atomic_fetch_and_add(volatile int64_t* _ptr, int64_t _addend) { return (int64_t)__sync_fetch_and_add(_ptr, _addend); }
inline int64_t atomic_compare_and_swap(volatile int64_t* _ptr, int64_t _old, int64_t _new) {
    return (int64_t)___sync_val_compare_and_swap(_ptr, _old, _new);
}

inline uint64_t atomic_fetch(volatile uint64_t* _ptr) { return (uint64_t)__sync_fetch_and_add(_ptr, 0); }
inline uint64_t atomic_fetch_and_set(volatile uint64_t* _ptr, uint64_t _val) {
    uint64_t result = (uint64_t)__sync_lock_test_and_set(_ptr, _val);
    __sync_lock_release(_ptr);
    return result;
}
inline uint64_t atomic_fetch_and_add(volatile uint64_t* _ptr, uint64_t _addend) { return (uint64_t)__sync_fetch_and_add(_ptr, _addend); }
inline uint64_t atomic_compare_and_swap(volatile uint64_t* _ptr, uint64_t _old, uint64_t _new) {
    return (uint64_t)___sync_val_compare_and_swap(_ptr, _old, _new);
}

#endif

inline int32_t atomic_test_and_add(volatile int32_t* _ptr, int32_t _test, int32_t _addend) {
    int32_t old_val;
    int32_t new_val = atomic_fetch(_ptr);
    do {
        old_val = new_val;
        new_val = atomic_compare_and_swap(_ptr, old_val, new_val >= _test ? _test : new_val + _addend);

    } while (old_val != new_val);

    return old_val;
}

inline uint32_t atomic_test_and_add(volatile uint32_t* _ptr, uint32_t _test, uint32_t _addend) {
    uint32_t old_val;
    uint32_t new_val = atomic_fetch(_ptr);
    do {
        old_val = new_val;
        new_val = atomic_compare_and_swap(_ptr, old_val, new_val >= _test ? _test : new_val + _addend);

    } while (old_val != new_val);

    return old_val;
}

inline int64_t atomic_test_and_add(volatile int64_t* _ptr, int64_t _test, int64_t _addend) {
    int64_t old_val;
    int64_t new_val = atomic_fetch(_ptr);
    do {
        old_val = new_val;
        new_val = atomic_compare_and_swap(_ptr, old_val, new_val >= _test ? _test : new_val + _addend);

    } while (old_val != new_val);

    return old_val;
}

inline uint64_t atomic_test_and_add(volatile uint64_t* _ptr, uint64_t _test, uint64_t _addend) {
    uint64_t old_val;
    uint64_t new_val = atomic_fetch(_ptr);
    do {
        old_val = new_val;
        new_val = atomic_compare_and_swap(_ptr, old_val, new_val >= _test ? _test : new_val + _addend);

    } while (old_val != new_val);

    return old_val;
}