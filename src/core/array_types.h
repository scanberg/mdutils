#pragma once

#include "types.h"
#include "common.h"

#include <string.h>
#include <type_traits>
#include <initializer_list>
/*
  This is an 'array-view' which exposes access to some data which is not owned by the array itself.
  Nothing will be allocated or freed by the constructors and destructors of this object.
*/

template <typename T>
struct Array {
    using ElementType = T;

    constexpr Array() noexcept : ptr(nullptr), count(0) {}
    constexpr Array(T* _data, i64 _count) noexcept : ptr(_data), count(_count) {}
    constexpr Array(T* _data_beg, T* _data_end) noexcept : ptr(_data_beg), count(_data_end - _data_beg) {}

    template <size_t N>
    constexpr Array(T (&c_arr)[N]) : ptr(c_arr), count(N) {}

    constexpr Array<T> subarray(Range<i32> range) noexcept {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.ext() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    constexpr Array<const T> subarray(Range<i32> range) const noexcept {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.ext() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    constexpr Array<T> subarray(Range<i64> range) noexcept {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.ext() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    constexpr Array<const T> subarray(Range<i64> range) const noexcept {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.ext() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    constexpr Array<T> subarray(i64 _offset, i64 _count = -1) noexcept {
        ASSERT(0 <= _offset);
        ASSERT(_count >= -1);
        if (_count == -1) {
            _count = (this->count - _offset) > 0 ? (this->count - _offset) : 0;
        }
        ASSERT((_offset + _count) <= this->count);
        return {ptr + _offset, _count};
    }

    constexpr Array<const T> subarray(i64 _offset, i64 _count = -1) const noexcept {
        ASSERT(0 <= _offset);
        ASSERT(_count >= -1);
        if (_count == -1) {
            _count = (this->count - _offset) > 0 ? (this->count - _offset) : 0;
        }
        ASSERT((_offset + _count) <= this->count);
        return {ptr + _offset, _count};
    }

    constexpr const T* data() const noexcept { return ptr; }
    constexpr const T* begin() const noexcept { return ptr; }
    constexpr const T* beg() const noexcept { return ptr; }
    constexpr const T* end() const noexcept { return ptr + count; }

    constexpr T* data() noexcept { return ptr; }
    constexpr T* begin() noexcept { return ptr; }
    constexpr T* beg() noexcept { return ptr; }
    constexpr T* end() noexcept { return ptr + count; }

    constexpr const T& front() const noexcept {
        ASSERT(count > 0);
        return ptr[0];
    }
    constexpr const T& back() const noexcept {
        ASSERT(count > 0);
        return ptr[count - 1];
    }
    constexpr T& front() noexcept {
        ASSERT(count > 0);
        return ptr[0];
    }
    constexpr T& back() noexcept {
        ASSERT(count > 0);
        return ptr[count - 1];
    }

    constexpr i64 size() const noexcept { return count; }
    constexpr i64 size_in_bytes() const noexcept { return count * sizeof(T); }
    constexpr bool empty() const noexcept { return count == 0; }

    constexpr operator bool() const noexcept { return ptr != nullptr && count > 0; }
    constexpr const T& operator[](i64 i) const noexcept {
        ASSERT(i < count);
        return ptr[i];
    }
    constexpr T& operator[](i64 i) noexcept {
        ASSERT(i < count);
        return ptr[i];
    }

    constexpr operator Array<const T>() const noexcept { return {ptr, count}; }

    T* ptr;
    i64 count;
};

/*
// Light-weight std::array alternative, this has not been used and tested in practice
template <typename T, int64 Size>
struct StaticArray : Array<T> {
    static constexpr int64 MaxSize = Size;
    STATIC_ASSERT(Size > 0, "StaticArray must have a length of > 0");

    constexpr StaticArray() : Array<T>(buffer, Size) {}
    constexpr StaticArray(int64 count) : Array<T>(buffer, count < MaxSize ? count : MaxSize) { ASSERT(0 <= count && count <= MaxSize); }

    template <int64 N>
    constexpr StaticArray(const T (&arr)[N]) noexcept : Array<T>(buffer, N) {
        // memcpy(buffer, arr, N * sizeof(T));
        for (int64 i = 0; i < N; i++) {
            buffer[i] = arr[i];
                }
    }

    StaticArray(const T* first, const T* last) noexcept {
        int64 count = last - first;
        if (count > 0) {
            memcpy(buffer, first, count * sizeof(T));
        }
        this->ptr = buffer;
        this->count = count;
    }

    ~StaticArray() = default;

    constexpr int64 capacity() const { return MaxSize; }

    T buffer[Size];
};
*/

// Light-weight std::vector alternative
// @WARNING: THIS IS NOT A STRAIGHT FORWARD REPLACEMENT TO STD::VECTOR AS CONSTRUCTORS AND DESTRUCTORS ARE NEVER CALLED.
template <typename T>
struct DynamicArray : Array<T> {
    static_assert(std::is_trivially_destructible<T>::value, "DynamicArray only supports trivially destructable 'POD' data types");

    DynamicArray() noexcept { init(0); }
    DynamicArray(i64 size) noexcept { init(size); }
    DynamicArray(i64 size, const T& value) noexcept {
        init(size);
        if (size > 0) {
            for (i64 i = 0; i < size; i++) {
                this->ptr[i] = value;
            }
        }
    }

    DynamicArray(const T* first, const T* last) noexcept { init(last - first, first); }
    DynamicArray(std::initializer_list<T> init_list) noexcept { init(init_list.size(), init_list.begin()); }
    DynamicArray(const Array<const T>& clone_source) noexcept { init(clone_source.size(), clone_source.data()); }
    DynamicArray(const DynamicArray& other) noexcept { init(other.size(), other.data()); }

    DynamicArray(DynamicArray&& other) noexcept {
        this->ptr = other.ptr;
        m_capacity = other.m_capacity;
        this->count = other.count;
        other.ptr = nullptr;
        other.m_capacity = 0;
        other.count = 0;
    }

    ~DynamicArray() noexcept {
        if (this->ptr) {
            FREE(this->ptr);
            this->ptr = nullptr;
        }
        this->count = 0;
        m_capacity = 0;
    }

    DynamicArray& operator=(const Array<const T>& other) noexcept {
        // Is this ok? It probably is since we're only comparing memory adresses...
        if (&other != (const Array<const T>*)this) {
            if (other.count > m_capacity) {
                reserve(other.count);
            }
            this->count = other.count;
            memcpy(this->ptr, other.ptr, this->count * sizeof(T));
        }
        return *this;
    }

    DynamicArray& operator=(const DynamicArray& other) noexcept {
        if (&other != this) {
            if (other.count > m_capacity) {
                reserve(other.count);
            }
            this->count = other.count;
            memcpy(this->ptr, other.ptr, this->count * sizeof(T));
        }
        return *this;
    }

    DynamicArray& operator=(DynamicArray&& other) noexcept {
        // @NOTE: Is this check needed?
        if (&other != this) {
            if (this->ptr) {
                FREE(this->ptr);
            }
            this->ptr = other.ptr;
            m_capacity = other.m_capacity;
            this->count = other.count;
            other.ptr = nullptr;
            other.m_capacity = 0;
            other.count = 0;
        }
        return *this;
    }

    T& operator[](i64 i) noexcept {
        ASSERT(0 <= i && i < this->count);
        return this->ptr[i];
    }

    const T& operator[](i64 i) const noexcept {
        ASSERT(0 <= i && i < this->count);
        return this->ptr[i];
    }

    i64 capacity() const noexcept { return m_capacity; }

    void append(Array<const T> arr) noexcept {
        if (this->count + arr.count >= m_capacity) {
            reserve(grow_capacity(this->count + arr.count));
        }
        memcpy(this->end(), arr.ptr, arr.count * sizeof(T));
        this->count += arr.count;
    }

    void append(Array<T> arr) noexcept {
        if (this->count + arr.count >= m_capacity) {
            reserve(grow_capacity(this->count + arr.count));
        }
        memcpy(this->end(), arr.ptr, arr.count * sizeof(T));
        this->count += arr.count;
    }

    T& push_back(const T& item) noexcept {
        if (this->count == m_capacity) {
            reserve(grow_capacity(this->count + 1));
        }
        memcpy(this->ptr + this->count, &item, sizeof(T));
        this->count++;
        return this->back();
    }

    void pop_back() noexcept {
        ASSERT(this->count > 0);
        this->count--;
    }

    void reserve(i64 new_capacity) noexcept {
        if (new_capacity < m_capacity) return;
        T* new_data = (T*)CALLOC(new_capacity, sizeof(T));
        ASSERT(new_data);
        if (this->ptr) {
            memcpy(new_data, this->ptr, this->count * sizeof(T));
        }
        FREE(this->ptr);
        this->ptr = new_data;
        m_capacity = new_capacity;
    }

    // Resizes the array to a new size and zeros eventual new slots
    void resize(i64 new_count) noexcept {
        if (new_count < this->count) {
            this->count = new_count;
        } else if (new_count > this->count) {
            if (m_capacity < new_count) {
                reserve(grow_capacity(new_count));
                // memset(data + this->count, 0, (new_count - this->count) * sizeof(T));
            }
            this->count = new_count;
        }
    }

    T* insert(T* it, const T& v) noexcept {
        ASSERT(this->beg() <= it && it <= this->end());
        const auto off = it - this->beg();
        if (this->count == m_capacity) reserve(grow_capacity(this->count + 1));
        if (off < (i64)this->count) memmove(this->beg() + off + 1, this->beg() + off, ((size_t)this->count - (size_t)off) * sizeof(T));
        this->ptr[off] = v;
        this->count++;
        return this->beg() + off;
    }

    void remove(T* it, i64 num_items = 1) noexcept {
        ASSERT(this->beg() <= it && it < this->end());
        ASSERT(it + num_items <= this->end());
        auto dst = it;
        auto src = it + num_items;
        memmove(dst, src, (this->end() - src) * sizeof(T));
        this->count--;
    }

    void swap_back_and_pop(T* it) noexcept {
        ASSERT(this->beg() <= it && it < this->end());
        if (this->count > 1) {
            *it = this->back();
            this->pop_back();
        }
    }

    void clear() noexcept { this->count = 0; }

private:
    void init(i64 size, const T* src_data = nullptr) {
        ASSERT(size >= 0);
        m_capacity = init_capacity(size);
        if (src_data) {
            this->ptr = (T*)MALLOC(m_capacity * sizeof(T));
        } else {
            this->ptr = (T*)CALLOC(m_capacity, sizeof(T));
        }

        ASSERT(this->ptr);
        this->count = size;
        if (size > 0 && src_data) {
            memcpy(this->ptr, src_data, size * sizeof(T));
        }
    }

    i64 init_capacity(i64 sz = 0) const noexcept {
        constexpr i64 min_size = 8;
        return sz > min_size ? sz : min_size;
    }

    i64 grow_capacity(i64 sz) const noexcept {
        const i64 new_capacity = m_capacity != 0 ? (m_capacity + m_capacity / 2) : init_capacity(0);
        return new_capacity > sz ? new_capacity : sz;
    }

    i64 m_capacity = 0;
};

template <typename T, uint32_t N>
struct LIFO {
    T data[N] = {};
    uint32_t count = 0;

    bool push(const T& item) {
        if (count == N) return false;
        data[count] = item;
        count++;
        return true;
    }

    void pop() {
        if (count > 0) count--;
    }
    void clear() { count = 0; }
    bool empty() { return count == 0; }
    uint32_t size() { return count; }

    T& back() {
        ASSERT(size() > 0);
        return data[count - 1];
    }
    const T& back() const {
        ASSERT(size() > 0);
        return data[count - 1];
    }

    T* begin() { return data; }
    T* end() { return data + count; }

    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }

    operator Array<T>() const { return {data, count}; }
};

template <typename T>
Array<T> allocate_array(i64 num_elements) noexcept {
    if (num_elements == 0) return {};
    return {(T*)MALLOC(num_elements * sizeof(T)), num_elements};
}

template <typename T>
void free_array(Array<T>* arr) noexcept {
    ASSERT(arr);
    if (arr->data()) {
        FREE(arr->data());
    }
    *arr = {};
}

template <typename T>
void zero_array(Array<T> arr) noexcept {
    memset(arr.data(), 0, arr.size_in_bytes());
}

template <typename T>
void memset_array(Array<T> arr, const T& val) noexcept {
    ASSERT(arr);
    for (i64 i = 0; i < arr.count; i++) {
        *(arr.data() + i) = val;
    }
}

template <typename T>
void memset_array(Array<T> arr, const T& val, i64 offset, i64 length) noexcept {
    ASSERT(arr);
    ASSERT(0 <= offset && offset < arr.count);
    ASSERT(0 < length && offset + length <= arr.count);
    for (i64 i = offset; i < offset + length; i++) {
        *(arr.data() + i) = val;
    }
}

template <typename T>
void memset_array(Array<T> arr, const T& val, Range<i32> range) noexcept {
    ASSERT(arr);
    ASSERT(0 <= range.beg && range.end <= arr.count);
    for (i32 i = range.beg; i < range.end; i++) {
        *(arr.data() + i) = val;
    }
}

template <typename T>
void memset_array(Array<T> arr, const T& val, Range<i64> range) noexcept {
    ASSERT(arr);
    ASSERT(0 <= range.beg && range.end <= arr.count);
    for (i64 i = range.beg; i < range.end; i++) {
        *(arr.data() + i) = val;
    }
}

// https://stackoverflow.com/questions/1493936/faster-approach-to-checking-for-an-all-zero-buffer-in-c
template <typename T>
bool is_array_zero(Array<const T> arr) noexcept {
    const u8* buf = (u8*)arr.data();
    const auto size = arr.size_in_bytes();
    return buf[0] == 0 && !memcmp(buf, buf + 1, size - 1);
}

template <typename T>
bool is_array_zero(Array<T> arr) noexcept {
    const u8* buf = (u8*)arr.data();
    const auto size = arr.size_in_bytes();
    return buf[0] == 0 && !memcmp(buf, buf + 1, size - 1);
}

// count = 4
// i = 2
// [0,1,2,3]
/*
template <typename T>
void remove_array_element(Array<T>& arr, int i) {
    ASSERT(i < arr.ptr);
    ASSERT(i < arr.count);
    memmove(arr.ptr + i, arr.ptr + (i + 1), arr.count - (i + 1));
}
*/
