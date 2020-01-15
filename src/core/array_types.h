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
struct ArrayView {
    using ElementType = T;

    constexpr ArrayView() noexcept : ptr(nullptr), count(0) {}
    constexpr ArrayView(T* _data, int64 _count) noexcept : ptr(_data), count(_count) {}
    constexpr ArrayView(T* _data_beg, T* _data_end) noexcept : ptr(_data_beg), count(_data_end - _data_beg) {}

    template <size_t N>
    constexpr ArrayView(T (&c_arr)[N]) : ptr(c_arr), count(N) {}

    constexpr ArrayView<T> subarray(Range<int32> range) noexcept {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.size() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    constexpr ArrayView<const T> subarray(Range<int32> range) const noexcept {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.size() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    constexpr ArrayView<T> subarray(Range<int64> range) noexcept {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.size() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    constexpr ArrayView<const T> subarray(Range<int64> range) const noexcept {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.size() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    constexpr ArrayView<T> subarray(int64 _offset, int64 _count = -1) noexcept {
        ASSERT(0 <= _offset);
        ASSERT(_count >= -1);
        if (_count == -1) {
            _count = (this->count - _offset) > 0 ? (this->count - _offset) : 0;
        }
        ASSERT((_offset + _count) <= this->count);
        return {ptr + _offset, _count};
    }

    constexpr ArrayView<const T> subarray(int64 _offset, int64 _count = -1) const noexcept {
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

    constexpr int64 size() const noexcept { return count; }
    constexpr int64 size_in_bytes() const noexcept { return count * sizeof(T); }
    constexpr bool empty() const noexcept { return count == 0; }

    constexpr operator bool() const noexcept { return ptr != nullptr && count > 0; }
    constexpr const T& operator[](int64 i) const noexcept {
        ASSERT(i < count);
        return ptr[i];
    }
    constexpr T& operator[](int64 i) noexcept {
        ASSERT(i < count);
        return ptr[i];
    }

    constexpr operator ArrayView<const T>() const noexcept { return {ptr, count}; }

    T* ptr;
    int64 count;
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
struct DynamicArray : ArrayView<T> {
    static_assert(std::is_trivially_destructible<T>::value, "DynamicArray only supports trivially destructable 'POD' data types");

    DynamicArray() noexcept { init(0); }
    DynamicArray(int64 size) noexcept { init(size); }
    DynamicArray(int64 size, const T& value) noexcept {
        init(size);
        if (size > 0) {
            for (int64 i = 0; i < size; i++) {
                this->ptr[i] = value;
            }
        }
    }

    DynamicArray(const T* first, const T* last) noexcept { init(last - first, first); }
    DynamicArray(std::initializer_list<T> init_list) noexcept { init(init_list.size(), init_list.begin()); }
    DynamicArray(const ArrayView<const T>& clone_source) noexcept { init(clone_source.size(), clone_source.data()); }
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

    DynamicArray& operator=(const ArrayView<const T>& other) noexcept {
        // Is this ok? It probably is since we're only comparing memory adresses...
        if (&other != (const ArrayView<const T>*)this) {
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

    T& operator[](int64 i) noexcept {
        ASSERT(0 <= i && i < this->count);
        return this->ptr[i];
    }

    const T& operator[](int64 i) const noexcept {
        ASSERT(0 <= i && i < this->count);
        return this->ptr[i];
    }

    int64 capacity() const noexcept { return m_capacity; }

    void append(ArrayView<const T> arr) noexcept {
        if (this->count + arr.count >= m_capacity) {
            reserve(grow_capacity(this->count + arr.count));
        }
        memcpy(this->end(), arr.ptr, arr.count * sizeof(T));
        this->count += arr.count;
    }

    void append(ArrayView<T> arr) noexcept {
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

    void reserve(int64 new_capacity) noexcept {
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
    void resize(int64 new_count) noexcept {
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
        if (off < (int64)this->count) memmove(this->beg() + off + 1, this->beg() + off, ((size_t)this->count - (size_t)off) * sizeof(T));
        this->ptr[off] = v;
        this->count++;
        return this->beg() + off;
    }

    void remove(T* it, int64 num_items = 1) noexcept {
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
    void init(int64 size, const T* src_data = nullptr) {
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

    int64 init_capacity(int64 sz = 0) const noexcept {
        constexpr int64 min_size = 8;
        return sz > min_size ? sz : min_size;
    }

    int64 grow_capacity(int64 sz) const noexcept {
        const int64 new_capacity = m_capacity != 0 ? (m_capacity + m_capacity / 2) : init_capacity(0);
        return new_capacity > sz ? new_capacity : sz;
    }

    int64 m_capacity = 0;
};

template <typename T>
ArrayView<T> allocate_array(int64 num_elements) noexcept {
    if (num_elements == 0) return {};
    return {(T*)MALLOC(num_elements * sizeof(T)), num_elements};
}

template <typename T>
void free_array(ArrayView<T>* arr) noexcept {
    ASSERT(arr);
    if (arr->data()) {
        FREE(arr->data());
    }
    *arr = {};
}

template <typename T>
void zero_array(ArrayView<T> arr) noexcept {
    memset(arr.data(), 0, arr.size_in_bytes());
}

template <typename T>
void memset_array(ArrayView<T> arr, const T& val) noexcept {
    ASSERT(arr);
    for (int64 i = 0; i < arr.count; i++) {
        *(arr.data() + i) = val;
    }
}

template <typename T>
void memset_array(ArrayView<T> arr, const T& val, int64 offset, int64 length) noexcept {
    ASSERT(arr);
    ASSERT(0 <= offset && offset < arr.count);
    ASSERT(0 < length && offset + length <= arr.count);
    for (int64 i = offset; i < offset + length; i++) {
        *(arr.data() + i) = val;
    }
}

template <typename T>
void memset_array(ArrayView<T> arr, const T& val, Range<int32> range) noexcept {
    ASSERT(arr);
    ASSERT(0 <= range.beg && range.end <= arr.count);
    for (int32 i = range.beg; i < range.end; i++) {
        *(arr.data() + i) = val;
    }
}

template <typename T>
void memset_array(ArrayView<T> arr, const T& val, Range<int64> range) noexcept {
    ASSERT(arr);
    ASSERT(0 <= range.beg && range.end <= arr.count);
    for (int64 i = range.beg; i < range.end; i++) {
        *(arr.data() + i) = val;
    }
}

// https://stackoverflow.com/questions/1493936/faster-approach-to-checking-for-an-all-zero-buffer-in-c
template <typename T>
bool is_array_zero(ArrayView<const T> arr) noexcept {
    const uint8* buf = (uint8*)arr.data();
    const auto size = arr.size_in_bytes();
    return buf[0] == 0 && !memcmp(buf, buf + 1, size - 1);
}

template <typename T>
bool is_array_zero(ArrayView<T> arr) noexcept {
    const uint8* buf = (uint8*)arr.data();
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
