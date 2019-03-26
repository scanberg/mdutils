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
	uint8 buffer[MaxSize] = {};
	// int32 length = 0;

	StringBuffer() {};

	template <int64 N>
	StringBuffer(const char(&cstr)[N]) {
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

	StringBuffer(StringBuffer&& other) {
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
			auto len = other.count < MaxSize ? other.count : MaxSize;
			memcpy(buffer, other.buffer, len);
			buffer[len - 1] = '\0';
		}
		return *this;
	}

	StringBuffer& operator=(const CString& cstr) {
		auto len = cstr.count < (MaxSize - 1) ? cstr.count : (MaxSize - 1);
		memcpy(buffer, cstr.beg(), len);
		buffer[len] = '\0';
		return *this;
	}

	StringBuffer& operator=(char c) {
		buffer[0] = c;
		buffer[1] = '\0';
		return *this;
	}

	template <int64 N>
	StringBuffer& operator=(const char(&cstr)[N]) {
		if (buffer != cstr) {
			constexpr auto len = N < MaxSize ? N : MaxSize;
			memcpy(buffer, cstr, len);
			buffer[len - 1] = '\0';
		}
		return *this;
	}

	StringBuffer& operator=(const char* cstr) {
		int64 len = (int64)strnlen(cstr, MaxSize);
		memcpy(buffer, cstr, len);
		buffer[len - 1] = '\0';
		return *this;
	}

	StringBuffer& operator+=(CString txt) {
		int64 offset = (int64)strnlen((const char*)buffer, MaxSize);
		int64 length = MaxSize - offset;
		if (length > 0) {
			strncpy((char*)buffer + offset, txt.cstr(), txt.count < length ? txt.count : length);
		}
		buffer[MaxSize - 1] = '\0';
		return *this;
	}

	char operator[](int64 i) const {
		ASSERT(i < MaxSize);
		return buffer[i];
	}

	char& operator[](int64 i) {
		ASSERT(i < MaxSize);
		return buffer[i];
	}

	//operator String() { return String(buffer, MaxSize); }
	//operator CString() const { return CString(buffer, strnlen((const char*)buffer, MaxSize)); }
	//operator char*() const { return (char*)buffer; }
	//operator const char*() const { return (const char*)buffer; }
	operator bool() const { return buffer[0] != '\0'; }

	int64 capacity() const { return MaxSize; }
	int64 size() const { return MaxSize; }

	char* cstr() const { return (char*)buffer; }

	const uint8* begin() const { return buffer; }
	const uint8* beg() const { return buffer; }
	const uint8* end() const { return buffer + MaxSize; }

	uint8* begin() { return buffer; }
	uint8* beg() { return buffer; }
	uint8* end() { return buffer + MaxSize; }
};


struct CString : Array<const uint8> {
    CString() = default;

	template <int64 N>
	CString(const StringBuffer<N>& buf) {
		ptr = buf.buffer;
		count = strnlen(buf.cstr(), buf.MaxSize);
	}

    CString(const char* cstr, int64 length = -1) {
        ptr = (const uint8*)cstr;
        if (length == -1) length = strlen(cstr);
        count = length;
    }

    CString(const uint8* cstr, int64 length = -1) {
        ptr = cstr;
        if (length == -1) length = strlen((const char*)cstr);
        count = length;
    }

    CString(const uint8* beg, const uint8* end) {
        ptr = beg;
        count = end - beg;
        ASSERT(count >= 0);
    }

    template <int64 length>
    CString(const char (&cstr)[length]) {
        ptr = (const uint8*)cstr;
        count = length;
    }

    CString substr(int64 _offset, int64 _count = -1) {
        auto array = subarray(_offset, _count);
        return {array.ptr, array.count};
    }

    const char* cstr() const { return (const char*)ptr; }
    int64 length() const { return count; }

    //operator const char*() { return (const char*)ptr; }
    operator bool() { return (ptr != 0 && count != 0); }
};

struct String : Array<uint8> {
    String() {
        ptr = 0;
        count = 0;
    }

	template <int64 N>
	String(StringBuffer<N>& buf) {
		ptr = buf.buffer;
		count = buf.MaxSize;
	}

    String(const String& other) : Array(other.ptr, other.count) {}

    String(char* cstr, int64 length) {
        ptr = (uint8*)cstr;
        count = length;
    }

    String(uint8* cstr, int64 length) {
        ptr = cstr;
        count = length;
    }

    String(uint8* beg, uint8* end) {
        ptr = beg;
        count = end - beg;
        ASSERT(count >= 0);
    }

    template <int64 length>
    String(char (&cstr)[length]) {
        ptr = (uint8*)cstr;
        count = length;
    }

    String substr(int64 _offset, int64 _count = -1) {
        auto array = subarray(_offset, _count);
        return {array.ptr, array.count};
    }

    char* cstr() const { return (char*)ptr; }
    int64 length() const { return count; }

    operator CString() { return CString(ptr, count); }
    //operator char*() { return (char*)ptr; }
    //operator const char*() { return (const char*)ptr; }
    operator bool() { return (ptr != 0 && count != 0); }
};

template <int64 N>
StringBuffer<N>::StringBuffer(const CString& str) {
	// @NOTE: MAX_LENGTH - 1 here because we copy from cstring which excludes \0
	auto len = str.length() < MaxSize - 1 ? str.length() : MaxSize - 1;
	memcpy(buffer, str.cstr(), len);
	buffer[len] = '\0';
}