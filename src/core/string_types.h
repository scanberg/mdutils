#pragma once

#include <string.h>

#include "types.h"
#include "array_types.h"

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

struct StringView;
struct CStringView;

namespace helper {
    constexpr i64 cexpr_strlen(const char* cstr) {
        i64 len = 0;
        while (*cstr != '\0') {
            ++cstr;
            ++len;
        }
        return len;
    }

    constexpr i64 cexpr_strnlen(const char* cstr, i64 max_length) {
        i64 len = 0;
        while (len < max_length && *cstr != '\0') {
            ++cstr;
            ++len;
        }
        return len;
    }

    constexpr i64 cexpr_strncpy(char* dst, const char* src, i64 max_length) {
        i64 i = 0;
        while (i < max_length && src[i] != '\0') {
            dst[i] = src[i];
            ++i;
        }
        return i;
    }

    template <typename T>
    constexpr void cexpr_copy(T* dst, const T* src, i64 count) {
        i64 i = 0;
        while (i < count) {
            dst[i] = src[i];
            i++;
        }
    }
}

// A buffer string for wrapping a char buffer[N]
template <i64 Size>
struct StringBuffer {
    static constexpr i64 MaxSize = Size;
    STATIC_ASSERT(MaxSize > 1, "Size of StringBuffer must be more than 1");
    char buffer[MaxSize] = {};

    constexpr StringBuffer() noexcept {};

    template <i64 N>
    constexpr StringBuffer(const char (&cstr)[N]) noexcept  {
        constexpr auto len = N < MaxSize ? N : MaxSize - 1;
        helper::cexpr_copy(buffer, cstr, len);
        //memcpy(buffer, cstr, len);
        buffer[len] = '\0';
    }

    constexpr StringBuffer(const char* cstr) noexcept {
        i64 len = helper::cexpr_strnlen(cstr, MaxSize - 1);
        helper::cexpr_copy(buffer, cstr, len);
        buffer[len] = '\0';
    }

    constexpr StringBuffer(const char* cstr, int64_t len) noexcept {
        ASSERT(len > 0);
        len = len < MaxSize - 1 ? len : MaxSize - 1;
        helper::cexpr_copy(buffer, cstr, len);
        buffer[len] = '\0';
    }

    constexpr StringBuffer(char c) noexcept {
        buffer[0] = c;
        buffer[1] = '\0';
    }

    constexpr StringBuffer(const StringBuffer& other) noexcept {
        helper::cexpr_copy(buffer, other.buffer, MaxSize);
        //memcpy(buffer, other.buffer, MaxSize);
        buffer[MaxSize - 1] = '\0';
    }

    template <i64 N>
    constexpr StringBuffer(const StringBuffer<N>& other) noexcept {
        constexpr auto len = N < MaxSize ? N : MaxSize;
        helper::cexpr_copy(buffer, other.buffer, len);
        buffer[len - 1] = '\0';
    }

    constexpr StringBuffer(StringBuffer&& other) noexcept {
        helper::cexpr_copy(buffer, other.buffer, MaxSize);
        buffer[MaxSize - 1] = '\0';
    }

    template <i64 N>
    constexpr StringBuffer(StringBuffer<N>&& other) noexcept {
        constexpr auto len = N < MaxSize ? N : MaxSize;
        helper::cexpr_copy(buffer, other.buffer, len);
        buffer[len - 1] = '\0';
    }

    StringBuffer(const CStringView& str);

    constexpr StringBuffer& operator=(const StringBuffer& other) noexcept {
        if (this != &other) {
            helper::cexpr_copy(buffer, other.buffer, MaxSize);
            buffer[MaxSize - 1] = '\0';
        }
        return *this;
    }

    template <i64 N>
    constexpr StringBuffer& operator=(const StringBuffer<N>& other) noexcept {
        if (this != &other) {
            constexpr auto len = N < (MaxSize - 1) ? N : (MaxSize - 1);
            helper::cexpr_copy(buffer, other.buffer, len);
            buffer[len] = '\0';
        }
        return *this;
    }

    StringBuffer& operator=(const CStringView& cstr);

    constexpr StringBuffer& operator=(char c) noexcept {
        buffer[0] = c;
        buffer[1] = '\0';
        return *this;
    }

    template <i64 N>
    constexpr StringBuffer& operator=(const char (&cstr)[N]) noexcept {
        if (buffer != cstr) {
            constexpr auto len = N < (MaxSize - 1) ? N : (MaxSize - 1);
            helper::cexpr_copy(buffer, cstr, len);
            buffer[len] = '\0';
        }
        return *this;
    }

    constexpr StringBuffer& operator=(const char* cstr) noexcept {
        const i64 len = helper::cexpr_strnlen(cstr, MaxSize - 1);
        helper::cexpr_copy(buffer, cstr, len);
        buffer[len] = '\0';
        return *this;
    }

    StringBuffer& operator+=(CStringView txt);

    constexpr char operator[](i64 i) const noexcept {
        ASSERT(0 <= i && i < MaxSize);
        return buffer[i];
    }

    constexpr char& operator[](i64 i) noexcept {
        ASSERT(0 <= i && i < MaxSize);
        return buffer[i];
    }

    // operator String() { return String(buffer, MaxSize); }
    // operator CString() const { return CString(buffer, strnlen((const char*)buffer, MaxSize)); }
    // operator char*() const { return (char*)buffer; }
    // operator const char*() const { return (const char*)buffer; }
    constexpr operator bool() const noexcept { return buffer[0] != '\0'; }

    constexpr i64 capacity() const noexcept { return MaxSize; }
    constexpr i64 length() const noexcept { return helper::cexpr_strnlen(buffer, MaxSize); }

    constexpr char* cstr() const noexcept { return (char*)buffer; }

    constexpr const char* begin() const noexcept { return buffer; }
    constexpr const char* beg() const noexcept { return buffer; }
    constexpr const char* end() const noexcept { return buffer + MaxSize; }

    constexpr char* begin() noexcept { return buffer; }
    constexpr char* beg() noexcept { return buffer; }
    constexpr char* end() noexcept { return buffer + MaxSize; }
};

struct CStringView : Array<const char> {
    constexpr CStringView() noexcept : Array() {};

    template <i64 N>
    constexpr CStringView(const StringBuffer<N>& buf) noexcept {
        ptr = buf.buffer;
        count = helper::cexpr_strnlen(buf.cstr(), buf.MaxSize);
    }

    template <i64 Length>
    constexpr CStringView(const char (&cstr)[Length]) noexcept {
        ptr = cstr;
        count = helper::cexpr_strnlen(cstr, Length - 1);
        ASSERT(count >= 0);
    }

    constexpr CStringView(const char* cstr, i64 len) noexcept {
        ptr = cstr;
        count = len;
        ASSERT(count >= 0);
    }

    constexpr CStringView(const char* beg, const char* end) noexcept {
        ptr = beg;
        count = end - beg;
        ASSERT(count >= 0);
    }

    constexpr CStringView(const char* cstr) noexcept {
        ptr = cstr;
        count = helper::cexpr_strlen(cstr);
    }

    constexpr CStringView substr(i64 _offset, i64 _count = -1) const noexcept {
        const auto arr = subarray(_offset, _count);
        return {arr.ptr, arr.count};
    }

    constexpr const char* cstr() const noexcept { return (const char*)ptr; }
    constexpr i64 length() const noexcept { return count; }

    //operator const char*() { return ptr; }
    //constexpr operator bool() noexcept { return (ptr != 0 && count != 0); }
};

struct StringView : Array<char> {
    constexpr StringView() = default;

    template <i64 N>
    constexpr StringView(StringBuffer<N>& buf) noexcept {
        ptr = buf.buffer;
        count = buf.MaxSize;
    }

    constexpr StringView(char* cstr, i64 length) noexcept {
        ptr = cstr;
        count = length;
    }

    constexpr StringView(char* beg, char* end) noexcept {
        ptr = beg;
        count = end - beg;
        ASSERT(count >= 0);
    }

    template <i64 length>
    constexpr StringView(char (&cstr)[length]) noexcept {
        ptr = cstr;
        count = length;
    }

    constexpr StringView substr(i64 _offset, i64 _count = -1) noexcept {
        auto arr = subarray(_offset, _count);
        return {arr.ptr, arr.count};
    }

    constexpr char* cstr() const noexcept { return ptr; }
    constexpr i64 length() const noexcept { return count; }

    constexpr operator CStringView() noexcept { return CStringView(ptr, count); }
    //operator char*() { return (char*)ptr; }
    //operator const char*() { return (const char*)ptr; }
    //constexpr operator bool() noexcept { return (ptr != 0 && count != 0); }
};

template <i64 N>
StringBuffer<N>::StringBuffer(const CStringView& str) {
    // @NOTE: MAX_LENGTH - 1 here because we copy from cstring which excludes \0
    auto len = str.length() < (MaxSize - 1) ? str.length() : (MaxSize - 1);
    helper::cexpr_copy(buffer, str.cstr(), len);
    buffer[len] = '\0';
}

template <i64 N>
StringBuffer<N>& StringBuffer<N>::operator=(const CStringView& str) {
    auto len = str.length() < (MaxSize - 1) ? str.length() : (MaxSize - 1);
    helper::cexpr_copy(buffer, str.cstr(), len);
    buffer[len] = '\0';
    return *this;
}

template <i64 N>
StringBuffer<N>& StringBuffer<N>::operator+=(CStringView txt) {
    i64 offset = (i64)helper::cexpr_strnlen((const char*)buffer, MaxSize);
    i64 length = MaxSize - offset;
    if (length > 0) {
        helper::cexpr_strncpy((char*)buffer + offset, txt.cstr(), txt.count < length ? txt.count : length);
    }
    buffer[MaxSize - 1] = '\0';
    return *this;
}

constexpr bool operator == (CStringView str_a, CStringView str_b) {
    if (str_a.size() != str_b.size()) return false;
    const char* ca = str_a.beg();
    const char* cb = str_b.beg();
    while (ca != str_a.end()) {
        if (*ca != *cb) return false;
        ++ca;
        ++cb;
    }
    return true;
}