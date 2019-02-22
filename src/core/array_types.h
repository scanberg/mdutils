#pragma once

#include "types.h"
#include "common.h"

#include <string.h>
#include <type_traits>

/*
  This is an 'array-view' which exposes access to some data which is not owned by the array itself.
  Nothing will be allocated or freed by the constructors and destructors of this object.
*/

template <typename T>
struct Array {
    Array() = default;
    Array(T* _data, int64 _count) : ptr(_data), count(_count) {}
    Array(T* _data_beg, T* _data_end) : ptr(_data_beg), count(_data_end - _data_beg) {}

    template <size_t N>
    Array(T (&c_arr)[N]) : ptr(c_arr), count(N) {}

    Array<T> subarray(int64 _offset, int64 _count = -1) {
        ASSERT(0 <= _offset);
        ASSERT(_count >= -1);
        if (_count == -1) {
            _count = (this->count - _offset) > 0 ? (this->count - _offset) : 0;
        }
        ASSERT((_offset + _count) <= this->count);
        return {ptr + _offset, _count};
    }

    Array<const T> subarray(int64 _offset, int64 _count = -1) const {
        ASSERT(0 <= _offset);
        ASSERT(_count >= -1);
        if (_count == -1) {
            _count = (this->count - _offset) > 0 ? (this->count - _offset) : 0;
        }
        ASSERT((_offset + _count) <= this->count);
        return {ptr + _offset, _count};
    }

    const T* data() const { return ptr; }
    const T* begin() const { return ptr; }
    const T* beg() const { return ptr; }
    const T* end() const { return ptr + count; }

    T* data() { return ptr; }
    T* begin() { return ptr; }
    T* beg() { return ptr; }
    T* end() { return ptr + count; }

    const T& front() const {
        ASSERT(count > 0);
        return ptr[0];
    }
    const T& back() const {
        ASSERT(count > 0);
        return ptr[count - 1];
    }
    T& front() {
        ASSERT(count > 0);
        return ptr[0];
    }
    T& back() {
        ASSERT(count > 0);
        return ptr[count - 1];
    }

    int64 size() const { return count; }
    int64 size_in_bytes() const { return count * sizeof(T); }

    operator bool() const { return ptr != nullptr && count > 0; }
    const T& operator[](int64 i) const {
        ASSERT(i < count);
        return ptr[i];
    }
    T& operator[](int64 i) {
        ASSERT(i < count);
        return ptr[i];
    }

    operator Array<const T>() const { return {ptr, count}; }

    T* ptr;
    int64 count;
};

// Light-weight std::array alternative
template <typename T, int64 Size>
struct StaticArray : Array<T> {
    static constexpr int64 MaxSize = Size;
    STATIC_ASSERT(Size > 0, "StaticArray must have a length of > 0");
    StaticArray() : Array<T>(buffer, Size) {}
    StaticArray(int64 count) : Array<T>(buffer, count < MaxSize ? count : MaxSize) { ASSERT(0 <= count && count <= MaxSize); }

    template <int64 N>
    StaticArray(const T (&arr)[N]) noexcept : Array<T>(buffer, N) {
        memcpy(buffer, arr, N * sizeof(T));
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

// Light-weight std::vector alternative
// @WARNING: THIS IS NOT A STRAIGHT FORWARD REPLACEMENT TO STD::VECTOR AS CONSTRUCTORS AND DESTRUCTORS ARE NEVER CALLED.
template <typename T>
struct DynamicArray : Array<T> {
    static_assert(std::is_trivially_destructible<T>::value, "DynamicArray only supports trivially destructable data types");

    static constexpr int64 INIT_CAPACITY = 8;
    DynamicArray() : m_capacity(INIT_CAPACITY) {
        this->ptr = (T*)CALLOC(m_capacity, sizeof(T));
        this->count = 0;
    }

    DynamicArray(int64 count) {
        count < 1 ? count = 1 : count;
        m_capacity = count;
        this->ptr = (T*)CALLOC(m_capacity, sizeof(T));
        this->count = m_capacity;
    }

    DynamicArray(int64 count, T value) : m_capacity(count > INIT_CAPACITY ? count : INIT_CAPACITY) {
        this->ptr = (T*)MALLOC(m_capacity * sizeof(T));
        this->count = count;
        for (int i = 0; i < count; i++) {
            this->ptr[i] = value;
        }
    }

    DynamicArray(const T* first, const T* last) noexcept {
        this->m_capacity = last - first;
        this->count = m_capacity;
        if (this->count > 0) {
            this->ptr = (T*)MALLOC(m_capacity * sizeof(T));
            memcpy(this->ptr, first, this->count * sizeof(T));
        }
    }

    DynamicArray(const Array<const T>& clone_source) : m_capacity(clone_source.count) {
        this->count = m_capacity;
        if (this->count > 0) {
            this->ptr = (T*)MALLOC(m_capacity * sizeof(T));
            memcpy(this->ptr, clone_source.ptr, this->count * sizeof(T));
        }
    }

    DynamicArray(const DynamicArray& other) : m_capacity(other.count) {
        this->count = m_capacity;
        if (this->count > 0) {
            this->ptr = (T*)MALLOC(m_capacity * sizeof(T));
            memcpy(this->ptr, other.ptr, this->count * sizeof(T));
        }
    }

    DynamicArray(DynamicArray&& other) {
        this->ptr = other.ptr;
        this->m_capacity = other.m_capacity;
        this->count = other.count;
        other.ptr = nullptr;
        other.m_capacity = 0;
        other.count = 0;
    }

    ~DynamicArray() {
        if (this->ptr) {
            FREE(this->ptr);
        }
        this->count = 0;
    }

    DynamicArray& operator=(const Array<const T>& other) {
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

    DynamicArray& operator=(const DynamicArray& other) {
        if (&other != this) {
            if (other.count > m_capacity) {
                reserve(other.count);
            }
            this->count = other.count;
            memcpy(this->ptr, other.ptr, this->count * sizeof(T));
        }
        return *this;
    }

    DynamicArray& operator=(DynamicArray&& other) {
        // @NOTE: Is this check needed?
        if (&other != this) {
            if (this->ptr) {
                FREE(this->ptr);
            }
            this->ptr = other.ptr;
            this->m_capacity = other.m_capacity;
            this->count = other.count;
            other.ptr = nullptr;
            other.m_capacity = 0;
            other.count = 0;
        }
        return *this;
    }

    int64 capacity() const { return m_capacity; }

    inline int64 _grow_capacity(int64 sz) const {
        int64 new_capacity = m_capacity ? (m_capacity + m_capacity / 2) : INIT_CAPACITY;
        return new_capacity > sz ? new_capacity : sz;
    }

    void append(Array<const T> arr) {
        if (this->count + arr.count >= m_capacity) {
            reserve(_grow_capacity(this->count + arr.count));
        }
        memcpy(this->end(), arr.ptr, arr.count * sizeof(T));
        this->count += arr.count;
    }

    void append(Array<T> arr) {
        if (this->count + arr.count >= m_capacity) {
            reserve(_grow_capacity(this->count + arr.count));
        }
        memcpy(this->end(), arr.ptr, arr.count * sizeof(T));
        this->count += arr.count;
    }

    T& push_back(const T& item) {
        if (this->count == m_capacity) {
            reserve(_grow_capacity(this->count + 1));
        }
        this->ptr[this->count] = item;
        this->count++;
        return this->back();
    }

    T pop_back() {
        ASSERT(this->count > 0);
        this->count--;
        return this->ptr[this->count];
    }

    void reserve(int64 new_capacity) {
        if (new_capacity < m_capacity) return;
        T* new_data = (T*)CALLOC(new_capacity, sizeof(T));
        if (this->ptr) {
            memcpy(new_data, this->ptr, this->count * sizeof(T));
        }
        FREE(this->ptr);
        this->ptr = new_data;
        m_capacity = new_capacity;
    }

    // Resizes the array to a new size and zeros eventual new slots
    void resize(int64 new_count) {
        if (new_count == this->count) {
            return;
        } else if (new_count < this->count) {
            this->count = new_count;
        } else {
            if (m_capacity < new_count) {
                reserve(_grow_capacity(new_count));
                // memset(this->data + this->count, 0, (new_count - this->count) * sizeof(T));
            }
            this->count = new_count;
        }
    }

    T* insert(T* it, const T& v) {
        ASSERT(this->beg() <= it && it <= this->end());
        const ptrdiff_t off = it - this->beg();
        if (this->count == m_capacity) reserve(_grow_capacity(this->count + 1));
        if (off < (int64)this->count) memmove(this->beg() + off + 1, this->beg() + off, ((size_t)this->count - (size_t)off) * sizeof(T));
        this->ptr[off] = v;
        this->count++;
        return this->beg() + off;
    }

    void remove(T* it, int64 num_items = 1) {
        ASSERT(this->beg() <= it && it < this->end());
        ASSERT(it + num_items <= this->end());
        auto dst = it;
        auto src = it + num_items;
        memmove(dst, src, (this->end() - src) * sizeof(T));
        this->count--;
    }

    void swap_back_and_pop(T* it) {
        ASSERT(this->beg() <= it && it < this->end());
        *it = this->back();
        pop_back();
    }

    void clear() { this->count = 0; }

private:
    int64 m_capacity;
};

template <typename T>
Array<T> allocate_array(int64 num_elements) {
    if (num_elements == 0) return {};
    return {(T*)MALLOC(num_elements * sizeof(T)), num_elements};
}

template <typename T>
void free_array(Array<T>* arr) {
    ASSERT(arr);
    if (arr->ptr) {
        FREE(arr->ptr);
    }
    arr->ptr = nullptr;
    arr->count = 0;
}

template <typename T>
void zero_array(Array<T> arr) {
    memset(arr.ptr, 0, arr.size_in_bytes());
}

template <typename T>
void memset_array(Array<T> arr, const T& val) {
    ASSERT(arr);
    for (int64 i = 0; i < arr.count; i++) {
        *(arr.ptr + i) = val;
    }
}

template <typename T>
void memset_array(Array<T> arr, const T& val, int64 offset, int64 length) {
    ASSERT(arr);
    ASSERT(0 <= offset && offset < arr.count);
    ASSERT(0 < length && offset + length <= arr.count);
    for (int64 i = offset; i < offset + length; i++) {
        *(arr.ptr + i) = val;
    }
}

// count = 4
// i = 2
// [0,1,2,3]
template <typename T>
void remove_array_element(Array<T>& arr, int i) {
    ASSERT(i < arr.ptr);
    ASSERT(i < arr.count);
    memmove(arr.ptr + i, arr.ptr + (i + 1), arr.count - (i + 1));
}
