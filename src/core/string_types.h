#pragma once

#include <string.h>

#include "types.h"
#include "array_types.h"

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

struct String;
struct CString;

// A buffer string for wrapping a char buffer[N]
template <int64 Size>
struct StringBuffer {
    static constexpr int64 MaxSize = Size;
    STATIC_ASSERT(MaxSize > 1, "Size of StringBuffer must be more than 1");
    char buffer[MaxSize] = {};
    // int32 length = 0;

    StringBuffer(){};

    template <int64 N>
    StringBuffer(const char (&cstr)[N]) {
        constexpr auto len = N < MaxSize ? N : MaxSize - 1;
        memcpy(buffer, cstr, len);
        buffer[len] = '\0';
    }

    StringBuffer(const char* cstr) {
        int64 len = (int64)strnlen(cstr, MaxSize - 1);
        memcpy(buffer, cstr, len);
        buffer[len] = '\0';
    }

    StringBuffer(char c) {
        buffer[0] = c;
        buffer[1] = '\0';
    }

    StringBuffer(const StringBuffer& other) {
        memcpy(buffer, other.buffer, MaxSize);
        buffer[MaxSize - 1] = '\0';
    }

    template <int64 N>
    StringBuffer(const StringBuffer<N>& other) {
        constexpr auto len = N < MaxSize ? N : MaxSize;
        memcpy(buffer, other.buffer, len);
        buffer[len - 1] = '\0';
    }

    StringBuffer(StringBuffer&& other) noexcept {
        memcpy(buffer, other.buffer, MaxSize);
        buffer[MaxSize - 1] = '\0';
    }

    template <int64 N>
    StringBuffer(StringBuffer<N>&& other) {
        constexpr auto len = N < MaxSize ? N : MaxSize;
        memcpy(buffer, other.buffer, len);
        buffer[len - 1] = '\0';
    }

    StringBuffer(const CString& str);

    StringBuffer& operator=(const StringBuffer& other) {
        if (this != &other) {
            memcpy(buffer, other.buffer, MaxSize);
            buffer[MaxSize - 1] = '\0';
        }
        return *this;
    }

    template <int64 N>
    StringBuffer& operator=(const StringBuffer<N>& other) {
        if (this != &other) {
            constexpr auto len = N < (MaxSize - 1) ? N : (MaxSize - 1);
            memcpy(buffer, other.buffer, len);
            buffer[len] = '\0';
        }
        return *this;
    }

    StringBuffer& operator=(const CString& cstr);

    StringBuffer& operator=(char c) {
        buffer[0] = c;
        buffer[1] = '\0';
        return *this;
    }

    template <int64 N>
    StringBuffer& operator=(const char (&cstr)[N]) {
        if (buffer != cstr) {
            constexpr auto len = N < (MaxSize - 1) ? N : (MaxSize - 1);
            memcpy(buffer, cstr, len);
            buffer[len] = '\0';
        }
        return *this;
    }

    StringBuffer& operator=(const char* cstr) {
        int64 len = (int64)strnlen(cstr, MaxSize);
        memcpy(buffer, cstr, len);
        buffer[len] = '\0';
        return *this;
    }

    StringBuffer& operator+=(CString txt);

    char operator[](int64 i) const {
        ASSERT(0 <= i && i < MaxSize);
        return buffer[i];
    }

    char& operator[](int64 i) {
        ASSERT(0 <= i && i < MaxSize);
        return buffer[i];
    }

    // operator String() { return String(buffer, MaxSize); }
    // operator CString() const { return CString(buffer, strnlen((const char*)buffer, MaxSize)); }
    // operator char*() const { return (char*)buffer; }
    // operator const char*() const { return (const char*)buffer; }
    operator bool() const { return buffer[0] != '\0'; }

    int64 capacity() const { return MaxSize; }
    int64 size() const { return MaxSize; }

    char* cstr() const { return (char*)buffer; }

    const char* begin() const { return buffer; }
    const char* beg() const { return buffer; }
    const char* end() const { return buffer + MaxSize; }

    char* begin() { return buffer; }
    char* beg() { return buffer; }
    char* end() { return buffer + MaxSize; }
};

struct CString : Array<const char> {
    constexpr CString() = default;

    template <int64 N>
    constexpr CString(const StringBuffer<N>& buf) noexcept {
        ptr = buf.buffer;
        count = strnlen(buf.cstr(), buf.MaxSize);
    }

    constexpr CString(const char* cstr, int64 len) noexcept {
        ptr = cstr;
        count = len;
        ASSERT(count >= 0);
    }

    constexpr CString(const char* beg, const char* end) noexcept {
        ptr = beg;
        count = end - beg;
        ASSERT(count >= 0);
    }

    template <int64 Length>
    constexpr CString(const char (&cstr)[Length]) noexcept {
        ptr = cstr;
        count = Length - 1;
        ASSERT(count >= 0);
    }

    constexpr CString substr(int64 _offset, int64 _count = -1) const noexcept {
        const auto arr = subarray(_offset, _count);
        return {arr.ptr, arr.count};
    }

    constexpr const char* cstr() const noexcept { return (const char*)ptr; }
    constexpr int64 length() const noexcept { return count; }

    operator const char*() { return ptr; }
    constexpr operator bool() noexcept { return (ptr != 0 && count != 0); }
};

struct String : Array<char> {
    constexpr String() = default;

    template <int64 N>
    constexpr String(StringBuffer<N>& buf) noexcept {
        ptr = buf.buffer;
        count = buf.MaxSize;
    }

    constexpr String(char* cstr, int64 length) noexcept {
        ptr = cstr;
        count = length;
    }

    constexpr String(char* beg, char* end) noexcept {
        ptr = beg;
        count = end - beg;
        ASSERT(count >= 0);
    }

    template <int64 length>
    constexpr String(char (&cstr)[length]) noexcept {
        ptr = cstr;
        count = length;
    }

    constexpr String substr(int64 _offset, int64 _count = -1) noexcept {
        auto arr = subarray(_offset, _count);
        return {arr.ptr, arr.count};
    }

    constexpr char* cstr() const noexcept { return ptr; }
    constexpr int64 length() const noexcept { return count; }

    constexpr operator CString() noexcept { return CString(ptr, count); }
    operator char*() { return (char*)ptr; }
    operator const char*() { return (const char*)ptr; }
    constexpr operator bool() noexcept { return (ptr != 0 && count != 0); }
};

template <int64 N>
StringBuffer<N>::StringBuffer(const CString& str) {
    // @NOTE: MAX_LENGTH - 1 here because we copy from cstring which excludes \0
    auto len = str.length() < (MaxSize - 1) ? str.length() : (MaxSize - 1);
    memcpy(buffer, str.cstr(), len);
    buffer[len] = '\0';
}

template <int64 N>
StringBuffer<N>& StringBuffer<N>::operator=(const CString& str) {
    auto len = str.length() < (MaxSize - 1) ? str.length() : (MaxSize - 1);
    memcpy(buffer, str.cstr(), len);
    buffer[len] = '\0';
    return *this;
}

template <int64 N>
StringBuffer<N>& StringBuffer<N>::operator+=(CString txt) {
    int64 offset = (int64)strnlen((const char*)buffer, MaxSize);
    int64 length = MaxSize - offset;
    if (length > 0) {
        strncpy((char*)buffer + offset, txt.cstr(), txt.count < length ? txt.count : length);
    }
    buffer[MaxSize - 1] = '\0';
    return *this;
}