#pragma once

#include <core/platform.h>
#include <atomic>

using std::atomic_int32_t;
using std::atomic_uint32_t;
using std::atomic_int64_t;
using std::atomic_uint64_t;

using std::atomic_load;
using std::atomic_store;
using std::atomic_fetch_add;
using std::atomic_fetch_sub;
using std::atomic_fetch_and;
using std::atomic_fetch_or;
using std::atomic_fetch_xor;
using std::atomic_exchange;

/*
#if PLATFORM_WINDOWS

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>

#define CACHE_LINE_SIZE 64

struct Semaphore {
	HANDLE id;
	atomic_uint32_t count;
};

namespace semaphore {
inline Semaphore create(uint32_t initial_count) {
    Semaphore sem;
    sem.id = CreateSemaphoreA(NULL, (LONG)initial_count, LONG_MAX, NULL);
    sem.count = initial_count;
}

inline Semaphore destroy(Semaphore* semaphore) {
	CloseHandle(semaphore->id);
}

inline bool release(Semaphore* semaphore, uint32_t ) {
	atomic_fetch_add(&semaphore->count, 1);
	if (ReleaseSemaphore(semaphore->id, 1, NULL) == FALSE) {
		atomic_fetch_sub(&semaphore->count, 1);
		return false;
	} else {
		return true;
	}
}

inline bool wait(Semaphore* semaphore, uint32_t milliseconds) {
	if (WaitForSingleObjectEx(semaphore->id, milliseconds, FALSE) == WAIT_OBJECT_0) {
		return true;
	} else {
		return false;
	}
}

inline bool try_aquire(Semaphore* semaphore) {
	return wait(semaphore, 0);
}

inline bool wait_aquire(Semaphore* semaphore) {
	return wait(semaphore, INFINITE);
}

inline uint32_t value(Semaphore* semaphore) {
	return atomic_load(&semaphore->count);
}


}

#endif
*/