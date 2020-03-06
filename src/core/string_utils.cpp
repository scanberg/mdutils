#include "string_utils.h"
#include <core/common.h>
#include <core/log.h>
#include <ctype.h>

#ifdef WIN32
#pragma warning(disable : 4996)  // Unsafe strncpy
#endif

#define MIN(x, y) ((x < y) ? (x) : (y))
#define MAX(x, y) ((x > y) ? (x) : (y))

// Extern memmem replacement, located in Railgun_Trolldom.c
extern "C" {
char* Railgun_Trolldom(char* pbTarget, char* pbPattern, uint32_t cbTarget, uint32_t cbPattern);
}

static inline bool internal_compare_ignore_case(const char* str_a, const char* str_b, i64 len) {
    for (i64 i = 0; i < len; i++) {
        if (tolower(str_a[i]) != tolower(str_b[i])) return false;
    }
    return true;
}

static inline bool internal_compare(const char* str_a, const char* str_b, i64 len) { return memcmp(str_a, str_b, len) == 0; }

bool compare(CStringView str_a, CStringView str_b) {
    if (str_a.count != str_b.count) return false;
    if (str_a.count == 0) return false;
    return internal_compare(str_a.ptr, str_b.ptr, str_a.count);
}

bool compare_ignore_case(CStringView str_a, CStringView str_b) {
    if (str_a.count != str_b.count) return false;
    if (str_a.count == 0) return false;
    return internal_compare_ignore_case(str_a.ptr, str_b.ptr, str_a.count);
}

bool compare_n(CStringView str_a, CStringView str_b, i64 num_chars) {
    const i64 len = MIN(str_a.count, str_b.count);
    // if (MIN(str_a.count, str_b.count) < num_chars) return false;
    return internal_compare(str_a.ptr, str_b.ptr, MIN(len, num_chars));
}

bool compare_n_ignore_case(CStringView str_a, CStringView str_b, i64 num_chars) {
    const i64 len = MIN(str_a.count, str_b.count);
    // if (len < num_chars) return false;
    return internal_compare_ignore_case(str_a.ptr, str_b.ptr, MIN(len, num_chars));
}

void copy(StringView dst, CStringView src) {
    ASSERT(dst.ptr);
    ASSERT(src.ptr);
    const auto len = MIN(dst.count, src.count);
    memcpy(dst.ptr, src.ptr, len);
    dst.ptr[src.count] = '\0';
}

void copy_n(StringView dst, CStringView src, i64 num_chars) {
    ASSERT(dst.ptr);
    ASSERT(src.ptr);
    const auto len = MIN(dst.count, MIN(src.count, num_chars));
    memcpy(dst.ptr, src.ptr, len);
    dst.ptr[num_chars] = '\0';
}

StringView allocate_string(CStringView str) {
    if (str.count == 0) return {};
    char* ptr = (char*)MALLOC(str.count + 1);
    if (!ptr) {
        LOG_ERROR("Could not allocate space for string.");
        return {};
    }
    strncpy((char*)ptr, str.cstr(), str.count + 1);
    return {ptr, str.count};
}

StringView allocate_string(i32 length) {
    if (length == 0) return {};
    char* ptr = (char*)MALLOC(length);
    if (!ptr) {
        LOG_ERROR("Could not allocate space for string.");
        return {};
    }
    return {ptr, length};
}

void free_string(StringView* str) {
    if (str->ptr) {
        FREE(str->ptr);
        str->ptr = nullptr;
        str->count = 0;
    }
}

CStringView peek_line(CStringView str) {
    if (str.length() == 0) {
        return {};
    }

    const char* line_beg = str.beg();
    const char* line_end = (const char*)memchr(str.beg(), '\n', str.length());
    if (!line_end) {
        line_end = str.end();
    }

    return {line_beg, line_end};
}

CStringView extract_line(CStringView& str) {
    const char* str_beg = str.ptr;
    const char* str_end = str.ptr + str.count;

    if (str_beg == str_end) {
        return {};
    }

    // If we start on a new line character or return character, skip these
    while (str_beg != str_end && (*str_beg == '\r' || *str_beg == '\n')) ++str_beg;

    const char* line_beg = str_beg;
    const char* line_end = (const char*)memchr(str_beg, '\n', str.length());
    if (!line_end) {
        line_end = str_end;
        str_beg = str_end;
    } else {
        line_end++;
        // Step over return and new line characters
        str_beg = MIN(line_end, str_end);
        while (str_beg != str_end && (*str_beg == '\r' || *str_beg == '\n')) ++str_beg;

        // Prune '/r and /n' from line
        while (line_end != line_beg && (*(line_end - 1) == '\r' || *(line_end - 1) == '\n')) --line_end;
    }

    str.ptr = str_beg;
    str.count = str_end - str_beg;

    return {line_beg, line_end};
}

/*
bool copy_line(String& line, CString& str) {
    const char* str_beg = str.data;
    const char* str_end = str.data + str.count;

    if (str_beg == str_end) {
        line = {};
        return false;
    }

    const char* line_beg = str_beg;
    const char* line_end = line_beg;

    // Find return or new line character
    while (line_end < str_end && *line_end != '\r' && *line_end != '\n') ++line_end;

    // Step over return or new line characters
    str_beg = MIN(line_end + 1, str_end);
    while (str_beg < str_end && (*str_beg == '\r' || *str_beg == '\n')) ++str_beg;

    // @NOTE: Do not modify line.count, its value contains the length of the buffer its pointing to
    auto count = MIN(line_end - line_beg, line.count - 1);
    line.data = (char*)memcpy(line.data, line_beg, count);
    line.data[count] = '\0';

    str.data = str_beg;
    str.count = str_end - str_beg;

    return true;
}
*/

ConversionResult<f32> to_float32(CStringView str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    f32 val = strtof(buf.cstr(), &end);
    return {val, end != buf.cstr()};
}

ConversionResult<f64> to_float64(CStringView str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    f64 val = strtod(buf.cstr(), &end);
    return {val, end != buf.cstr()};
}

ConversionResult<i32> to_int32(CStringView str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    i32 val = strtol(buf.cstr(), &end, 10);
    return {val, end != buf.cstr()};
}

ConversionResult<i64> to_int64(CStringView str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    i64 val = strtoll(buf.cstr(), &end, 10);
    return {val, end != buf.cstr()};
}

CStringView trim(CStringView str) {
    const char* beg = str.ptr;
    const char* end = str.ptr + str.count;

    while (beg < end && is_whitespace(*beg)) ++beg;
    while (end > beg && (is_whitespace(*(end - 1)) || *(end - 1) == '\0')) --end;

    return CStringView(beg, end - beg);
}

StringView trim(StringView str) {
    char* beg = str.ptr;
    char* end = str.ptr + str.count;

    while (beg < end && is_whitespace(*beg)) ++beg;
    while (end > beg && is_whitespace(*(end - 1))) --end;

    return StringView(beg, end - beg);
}

StringView allocate_and_read_textfile(CStringView filename) {
    FILE* file = fopen(filename, "rb");
    defer {
        if (file) fclose(file);
    };

    if (!file) return {};

    fseeki64(file, 0, SEEK_END);
    i64 file_size = ftelli64(file);
    rewind(file);

    if (file_size <= 0) return {};

    char* ptr = (char*)MALLOC(file_size + 1);
    if (!ptr) {
        LOG_ERROR("Could not allocate space for string.");
        return {};
    }
    fread(ptr, 1, file_size, file);
    ptr[file_size] = '\0';

    return {ptr, file_size + 1};
}

i64 find_pattern_in_file(CStringView filename, CStringView pattern) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        LOG_ERROR("Could not open file '.*s'", filename.length(), filename.beg());
        return -1;
    }
    defer { fclose(file); };

    constexpr i64 buf_size = MEGABYTES(1);
    char* buf = (char*)TMP_MALLOC(buf_size);
    defer { TMP_FREE(buf); };

    i64 byte_offset = 0;
    i64 bytes_read = (i64)fread(buf, 1, buf_size, file);

    CStringView str = {buf, (i64)bytes_read};
    if (CStringView match = find_string(str, pattern)) {
        return (i64)(match.beg() - buf);
    }

    const i64 chunk_size = buf_size - pattern.size_in_bytes();
    while (bytes_read == buf_size) {
        // Copy potential 'cut' pattern at end of buffer
        memcpy(buf, buf + buf_size - pattern.size_in_bytes(), pattern.size_in_bytes());
        bytes_read = (i64)fread(buf + pattern.size_in_bytes(), 1, chunk_size, file) + pattern.size_in_bytes();
        byte_offset += chunk_size;

        str = {buf, (i64)bytes_read};
        if (CStringView match = find_string(str, pattern)) {
            return byte_offset + (i64)(match.beg() - buf);
        }
    }
    return -1;
}

// Finds all occurrences with offsets (in bytes) of a pattern within a file
DynamicArray<i64> find_patterns_in_file(CStringView filename, CStringView pattern) {
    DynamicArray<i64> offsets;

    FILE* file = fopen(filename, "rb");
    if (!file) {
        LOG_ERROR("Could not open file '.*s'", filename.length(), filename.beg());
        return offsets;
    }
    defer { fclose(file); };

    constexpr i64 buf_size = MEGABYTES(1);
    char* buf = (char*)TMP_MALLOC(buf_size);
    defer { TMP_FREE(buf); };

    i64 byte_offset = 0;
    i64 bytes_read = (i64)fread(buf, 1, buf_size, file);

    CStringView str = {buf, (i64)bytes_read};
    while (CStringView match = find_string(str, pattern)) {
        offsets.push_back((i64)(match.beg() - buf));
        str = {match.end(), str.end()};
    }

    const i64 chunk_size = buf_size - pattern.size_in_bytes();
    while (bytes_read == buf_size) {
        // Copy potential 'cut' pattern at end of buffer
        memcpy(buf, buf + buf_size - pattern.size_in_bytes(), pattern.size_in_bytes());
        bytes_read = (i64)fread(buf + pattern.size_in_bytes(), 1, chunk_size, file) + pattern.size_in_bytes();
        byte_offset += chunk_size;

        str = {buf, (i64)bytes_read};
        while (CStringView match = find_string(str, pattern)) {
            offsets.push_back(byte_offset + (i64)(match.beg() - buf));
            str = {match.end(), str.end()};
        }
    }

    return offsets;
}

CStringView get_directory(CStringView url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);

    const char* beg = url.begin();
    const char* end = url.end() - 1;

    while (end != beg && *end != '\\' && *end != '/') {
        end--;
    }

    return CStringView(beg, end - beg);
}

CStringView get_file(CStringView url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);

    const char* beg = url.end() - 1;
    const char* end = url.end();

    while (beg != url.begin() && *beg != '\\' && *beg != '/') {
        beg--;
    }
    if (*beg == '\\' || *beg == '/') beg++;

    return CStringView(beg, end - beg);
}

CStringView get_file_without_extension(CStringView url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);

    const char* beg = url.end() - 1;
    const char* end = url.end();

    while (beg != url.begin() && *beg != '\\' && *beg != '/') beg--;
    if (beg != url.begin()) beg++;

    while (end != beg && *end != '.') end--;

    return CStringView(beg, end - beg);
}

CStringView get_file_extension(CStringView url) {
    if (url.count == 0) {
        return {};
    }

    url = trim(url);

    const char* beg = url.end() - 1;
    const char* end = url.end();

    while (beg != url.begin() && *beg != '.' && *beg != '\\' && *beg != '/') beg--;

    if (beg == url.begin() || *beg == '\\' || *beg == '/') {
        return {};
    }

    beg++;  // skip '.'
    return {beg, end - beg};
}

StringView get_file_extension(StringView url) {
    if (url.count == 0) {
        return {};
    }

    url = trim(url);

    char* beg = url.end() - 1;
    char* end = url.end();

    while (beg != url.begin() && *beg != '.' && *beg != '\\' && *beg != '/') beg--;

    if (beg == url.begin() || *beg == '\\' || *beg == '/') {
        return {};
    }

    beg++;  // skip '.'
    return {beg, end - beg};
}

inline static bool char_in_string(char c, CStringView str) {
    for (i64 i = 0; i < str.count; i++) {
        if (c == str[i]) return true;
    }
    return false;
}

StringBuffer<256> get_relative_path(CStringView from, CStringView to) {
    const char* c_from = from.beg();
    const char* c_to = to.beg();
    while (c_from != from.end() && c_to != to.end() && *c_from == *c_to) {
        c_from++;
        c_to++;
    }

    // If they have nothing in common. Return absolute path of to
    if (c_to == to.beg()) {
        return to;
    }

    int dir_count = 0;
    for (const char* c = c_from; c != from.end(); c++ /* <- LOL! */) {
        if (*c == '\\' || *c == '/') dir_count++;
    }

    StringBuffer<256> res;
    int offset = 0;
    for (int i = 0; i < dir_count; i++) {
        offset += snprintf(res.cstr() + offset, res.capacity(), "../");
    }

    StringBuffer<256> to_buf = CStringView(c_to, to.end());
    snprintf(res.cstr() + offset, res.capacity(), "%s", to_buf.beg());

    return res;
}

StringBuffer<256> get_absolute_path(CStringView absolute_reference, CStringView relative_file) {
    CStringView abs_dir = get_directory(absolute_reference);
    if (relative_file.count < 3) return {};

    StringBuffer<256> res;
    // If relative path is really an absolute path, just return that
    if (relative_file[0] == '/' || relative_file[1] == ':') {
        res = relative_file;
        return res;
    }

    int dir_count = 0;
    for (const char* c = relative_file.beg(); c < relative_file.end(); c += 3) {
        if (c[0] == '.' && (c + 1) != relative_file.end() && c[1] == '.' && (c + 2) != relative_file.end() && (c[2] == '/' || c[2] == '\\'))
            dir_count++;
    }

    const char* c = abs_dir.end() - 1;
    while (c > abs_dir.beg() && dir_count > 0) {
        if (*c == '/' || *c == '\\') {
            if (--dir_count == 0) break;
        }
        c--;
    }
    if (dir_count > 0 || c == abs_dir.beg()) return res;

    CStringView base_dir(abs_dir.beg(), c + 1);
    res = base_dir;
    StringBuffer<128> file = relative_file;
    snprintf(res.cstr() + base_dir.size(), res.capacity() - base_dir.size(), "/%s", file.beg());

    return res;
}

void convert_backslashes(StringView str) {
    for (char* c = str.beg(); c != str.end(); c++) {
        if (*c == '\\') *c = '/';
    }
}

CStringView extract_parentheses(CStringView str) {
    const char* beg = str.beg();

    while (beg != str.end() && *beg != '(') beg++;
    if (beg == str.end()) return {beg, str.end()};

    const char* end = beg + 1;
    int count = 1;
    while (end != str.end()) {
        if (*end == '(')
            count++;
        else if (*end == ')' && --count == 0) {
            end++;
            break;
        }
        end++;
    }

    return {beg, end};
}

CStringView extract_parentheses_contents(CStringView str) {
    CStringView p = extract_parentheses(str);
    if (p.count < 2) return p;
    return {p.beg() + 1, p.end() - 1};
}

CStringView find_string(CStringView target, CStringView pattern) {
    if (pattern.count == 0 || target.count < pattern.count) return {};

    char* ptr = Railgun_Trolldom((char*)target.cstr(), (char*)pattern.cstr(), (u32)target.size_in_bytes(), (u32)pattern.size_in_bytes());
    if (ptr) {
        return {(char*)ptr, pattern.length()};
    }
    return {};
}

DynamicArray<CStringView> tokenize(CStringView str, char delimiter) {
    DynamicArray<CStringView> tokens;

    const char* beg = str.beg();
    const char* end = str.beg();

    while (end != str.end() && *end != '\0') {
        while (end != str.end() && *end != '\0' && *end != delimiter) end++;
        tokens.push_back(CStringView(beg, end));
        beg = end;
        while (beg != str.end() && *beg != '\0' && *beg == delimiter) beg++;
        end = beg;
    }

    return tokens;
}

DynamicArray<CStringView> tokenize(CStringView str, CStringView delimiters) {
    DynamicArray<CStringView> tokens;

    const char* beg = str.beg();
    const char* end = str.beg();

    while (end != str.end() && *end != '\0') {
        while (end != str.end() && *end != '\0' && !char_in_string(*end, delimiters)) end++;
        tokens.push_back(CStringView(beg, end));
        beg = end;
        while (beg != str.end() && *end != '\0' && char_in_string(*beg, delimiters)) beg++;
        end = beg;
    }

    return tokens;
}

constexpr char delimiter = ':';
constexpr char wildcard = '*';

// Range extraction functionality
bool is_range(CStringView arg) {
    for (const char* c = arg.beg(); c != arg.end(); c++) {
        if (is_digit(*c)) continue;
        if (*c == delimiter) return true;
        if (*c == wildcard) return true;
    }
    return false;
}

bool extract_range(Range<i32>* range, CStringView arg) {
    if (arg.count == 0) {
        *range = {-1, -1};
        return false;
    }

    if (arg.count == 1 && arg[0] == wildcard) {
        range->beg = -1;
        range->end = -1;
        return true;
    }

    const char* mid = arg.beg();
    while (mid != arg.end() && *mid != delimiter) mid++;
    if (mid == arg.end()) return false;

    CStringView str_first(arg.beg(), mid);
    CStringView str_last(mid + 1, arg.end());

    if (str_first.count == 1 && str_first[0] == wildcard) {
        range->beg = -1;
    } else {
        auto res = to_int32(str_first);
        if (!res) return false;
        range->beg = res;
    }

    if (str_last.count == 1 && str_last[0] == wildcard) {
        range->end = -1;
    } else {
        auto res = to_int32(str_last);
        if (!res) return false;
        range->end = res;
    }

    return true;
}

bool extract_ranges(DynamicArray<Range<i32>>* ranges, Array<const CStringView> args) {
    ASSERT(ranges);

    for (auto arg : args) {
        if (is_range(arg)) {
            Range<i32> r;
            if (!extract_range(&r, arg)) return false;
            ranges->push_back(r);
        } else {
            auto res = to_int(arg);
            if (!res.success) return false;
            ranges->push_back({res.value, res.value});
        }
    }
    return true;
}

void print_string(CStringView str) { printf("%.*s", (int)str.count, str.ptr); }