#pragma once

#include <stdint.h>

using i8 = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float f32;
typedef double f64;

template <typename T>
struct Range {
    union {
        T beg = 0, min;
    };
    union {
        T end = 0, max;
    };

    constexpr Range() = default;
    constexpr Range(T lo, T hi) : beg(lo), end(hi){};

    template <typename U>
    constexpr Range(const Range<U>& other)
        : beg((T)other.beg), end((T)other.end) {}

    constexpr operator bool() const { return beg != end && beg < end; }
	constexpr Range& operator +=(T val) {
		beg += val;
		end += val;
		return *this;
	}
	constexpr Range& operator -=(T val) {
		beg -= val;
		end -= val;
		return *this;
	}
    constexpr T ext() const { return end - beg; }
};

template <typename T>
constexpr bool operator==(const Range<T>& r_a, const Range<T>& r_b) {
    return r_a.beg == r_b.beg && r_a.end == r_b.end;
}

template <typename T>
constexpr bool operator!=(const Range<T>& r_a, const Range<T>& r_b) {
    return !(r_a == r_b);
}

template <typename T>
constexpr Range<T> operator+(const Range<T>& range, T val) {
	return { range.x + val, range.y + val };
}

template <typename T>
constexpr Range<T> operator-(const Range<T>& range, T val) {
	return { range.x - val, range.y - val };
}

template <typename T>
constexpr bool ranges_overlap(const Range<T>& a, const Range<T>& b) {
    return (a.beg < b.end && b.beg < a.end);
}