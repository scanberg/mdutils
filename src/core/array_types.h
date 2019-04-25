#pragma once

#include "types.h"
#include "common.h"

#include <string.h>
#include <type_traits>

template <class T>
struct InitializerList {
    using value_type = T;
    using reference = const T&;
    using const_reference = const T&;

    using iterator = const T*;
    using const_iterator = const T*;

    constexpr InitializerList() noexcept = default;
    constexpr InitializerList(const T* first_arg, const T* last_arg) noexcept : first(first_arg), last(last_arg) {}

    constexpr const T* beg() const noexcept { return first; }
    constexpr const T* begin() const noexcept { return first; }
    constexpr const T* end() const noexcept { return last; }
    constexpr int64 size() const noexcept { return (int64)(_Last - _First); }

    const T* first = nullptr;
    const T* last = nullptr;
};

template <class T>
constexpr const T* begin(InitializerList<T> list) noexcept {
    return list.begin();
}

template <class T>
constexpr const T* end(InitializerList<T> list) noexcept {
    return list.end();
}

/*
  This is an 'array-view' which exposes access to some data which is not owned by the array itself.
  Nothing will be allocated or freed by the constructors and destructors of this object.
*/

template <typename T>
struct Array {
    using ElementType = T;

    Array() = default;
    Array(T* _data, int64 _count) : ptr(_data), count(_count) {}
    Array(T* _data_beg, T* _data_end) : ptr(_data_beg), count(_data_end - _data_beg) {}

    template <size_t N>
    Array(T (&c_arr)[N]) : ptr(c_arr), count(N) {}

    Array<T> subarray(Range<int32> range) {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.size() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    Array<const T> subarray(Range<int32> range) const {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.size() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    Array<T> subarray(Range<int64> range) {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.size() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

    Array<const T> subarray(Range<int64> range) const {
        ASSERT(0 <= range.beg);
        ASSERT(range.end <= this->count);
        ASSERT(range.size() >= 0);
        return {ptr + range.beg, range.end - range.beg};
    }

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

    T* ptr = nullptr;
    int64 count = 0;
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
    static_assert(std::is_trivially_destructible<T>::value, "DynamicArray only supports trivially destructable 'POD' data types");

    DynamicArray() noexcept : m_capacity(init_capacity()) {
        this->ptr = (T*)CALLOC(m_capacity, sizeof(T));
        this->count = 0;
    }

    DynamicArray(int64 size) noexcept {
		init(size);
    }

    DynamicArray(int64 size, const T& value) noexcept {
        init(size);
        if (size > 0) {
            for (int64 i = 0; i < size; i++) {
                this->ptr[i] = value;
            }
        }
    }

    DynamicArray(const T* first, const T* last) noexcept {
        const auto size = last - first;
        init(size, first);
    }

    DynamicArray(InitializerList<T> init_list) noexcept {
        init(init_list.size(), init_list.begin());
    }

    DynamicArray(const Array<const T>& clone_source) noexcept {
		init(clone_source.size(), clone_source.data());
    }

    DynamicArray(const DynamicArray& other) noexcept {
		init(other.size(), other.data());
    }

    DynamicArray(DynamicArray&& other) noexcept {
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
            this->ptr = nullptr;
        }
        this->count = 0;
        this->m_capacity = 0;
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
            this->m_capacity = other.m_capacity;
            this->count = other.count;
            other.ptr = nullptr;
            other.m_capacity = 0;
            other.count = 0;
        }
        return *this;
    }

    int64 capacity() const noexcept { return m_capacity; }

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
        this->ptr[this->count] = item;
        this->count++;
        return this->back();
    }

    T pop_back() noexcept {
        ASSERT(this->count > 0);
        this->count--;
        return this->ptr[this->count];
    }

    void reserve(int64 new_capacity) noexcept {
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
    void resize(int64 new_count) noexcept {
        if (new_count == this->count) {
            return;
        } else if (new_count < this->count) {
            this->count = new_count;
        } else {
            if (m_capacity < new_count) {
                reserve(grow_capacity(new_count));
                // memset(this->data + this->count, 0, (new_count - this->count) * sizeof(T));
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
        if (it != &this->back()) *it = this->back();
        pop_back();
    }

    void clear() noexcept { this->count = 0; }

private:
	inline void init(int64 size, T* src_data = nullptr) {
        ASSERT(size >= 0);
        this->m_capacity = init_capacity(size);
        if (size > 0) {
            this->ptr = (T*)MALLOC(this->m_capacity * sizeof(T));
            ASSERT(this->ptr);
            if (src_data) {
                memcpy(this->ptr, src_data, size * sizeof(T));
            }
        }
        this->count = size;
	}

    inline int64 init_capacity(int64 sz = 0) const noexcept {
        constexpr int64 min_size = 8;
        return sz > min_size ? sz : min_size;
    }

    inline int64 grow_capacity(int64 sz) const noexcept {
        const int64 new_capacity = m_capacity != 0 ? (m_capacity + m_capacity / 2) : init_capacity(0);
        return new_capacity > sz ? new_capacity : sz;
    }

    int64 m_capacity = 0;
};

template <typename T>
Array<T> allocate_array(int64 num_elements) noexcept {
    if (num_elements == 0) return {};
    return {(T*)MALLOC(num_elements * sizeof(T)), num_elements};
}

template <typename T>
void free_array(Array<T>* arr) noexcept {
    ASSERT(arr);
    if (arr->ptr) {
        FREE(arr->ptr);
    }
    arr->ptr = nullptr;
    arr->count = 0;
}

template <typename T>
void zero_array(Array<T> arr) noexcept {
    memset(arr.ptr, 0, arr.size_in_bytes());
}

template <typename T>
void memset_array(Array<T> arr, const T& val) noexcept {
    ASSERT(arr);
    for (int64 i = 0; i < arr.count; i++) {
        *(arr.ptr + i) = val;
    }
}

template <typename T>
void memset_array(Array<T> arr, const T& val, int64 offset, int64 length) noexcept {
    ASSERT(arr);
    ASSERT(0 <= offset && offset < arr.count);
    ASSERT(0 < length && offset + length <= arr.count);
    for (int64 i = offset; i < offset + length; i++) {
        *(arr.ptr + i) = val;
    }
}

template <typename T>
void memset_array(Array<T> arr, const T& val, Range<int32> range) noexcept {
    ASSERT(arr);
    ASSERT(0 <= range.beg && range.end <= arr.count);
    for (int32 i = range.beg; i < range.end; i++) {
        *(arr.ptr + i) = val;
    }
}

template <typename T>
void memset_array(Array<T> arr, const T& val, Range<int64> range) noexcept {
    ASSERT(arr);
    ASSERT(0 <= range.beg && range.end <= arr.count);
    for (int64 i = range.beg; i < range.end; i++) {
        *(arr.ptr + i) = val;
    }
}

// https://stackoverflow.com/questions/1493936/faster-approach-to-checking-for-an-all-zero-buffer-in-c
template <typename T>
bool is_array_zero(Array<const T> arr) noexcept {
    const uint8* buf = (uint8*)arr.data();
    const auto size = arr.size_in_bytes();
    return buf[0] == 0 && !memcmp(buf, buf + 1, size - 1);
}

template <typename T>
bool is_array_zero(Array<T> arr) noexcept {
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
