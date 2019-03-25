#include <core/string_utils.h>
#include <core/common.h>
#include <stdio.h>
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

static inline bool internal_compare_ignore_case(const uint8* str_a, const uint8* str_b, int64 len) {
    for (int64 i = 0; i < len; i++) {
        if (tolower(str_a[i]) != tolower(str_b[i])) return false;
    }
    return true;
}

static inline bool internal_compare(const uint8* str_a, const uint8* str_b, int64 len) { return memcmp(str_a, str_b, len) == 0; }

bool compare(CString str_a, CString str_b) {
    if (str_a.count != str_b.count) return false;
    if (str_a.count == 0) return false;
    return internal_compare(str_a.ptr, str_b.ptr, str_a.count);
}

bool compare_ignore_case(CString str_a, CString str_b) {
    if (str_a.count != str_b.count) return false;
    if (str_a.count == 0) return false;
    return internal_compare_ignore_case(str_a.ptr, str_b.ptr, str_a.count);
}

bool compare_n(CString str_a, CString str_b, int64 num_chars) {
    const int64 len = MIN(str_a.count, MIN(str_b.count, num_chars));
    if (len < num_chars) return false;
    return internal_compare(str_a.ptr, str_b.ptr, len);
}

bool compare_n_ignore_case(CString str_a, CString str_b, int64 num_chars) {
    const int64 len = MIN(str_a.count, MIN(str_b.count, num_chars));
    if (len < num_chars) return false;
    return internal_compare_ignore_case(str_a.ptr, str_b.ptr, len);
}

void copy(String dst, CString src) {
    ASSERT(dst.ptr != 0);
    ASSERT(src.ptr != 0);
    auto len = MIN(dst.count, src.count);
    memcpy(dst.ptr, src.ptr, len);
    dst.ptr[src.count] = '\0';
}

void copy_n(String dst, CString src, int64 num_chars) {
    ASSERT(dst.ptr != 0);
    ASSERT(src.ptr != 0);
    auto len = MIN(dst.count, src.count);
    len = MIN(len, num_chars);
    memcpy(dst.ptr, src.ptr, len);
    dst.ptr[num_chars] = '\0';
}

String allocate_string(CString str) {
    if (str.count == 0) return {};
    uint8* ptr = (uint8*)MALLOC(str.count + 1);
    strncpy((char*)ptr, str.cstr(), str.count + 1);
    return {ptr, str.count};
}

String allocate_string(int32 length) {
    if (length == 0) return {};
    uint8* ptr = (uint8*)MALLOC(length);
    return {ptr, length};
}

void free_string(String* str) {
    if (str->ptr) {
        FREE(str->ptr);
        str->ptr = nullptr;
        str->count = 0;
    }
}

CString peek_line(CString str) {
    if (str.length() == 0) {
        return {};
    }

    const uint8* line_beg = str.beg();
    const uint8* line_end = (const uint8*)memchr(str.beg(), '\n', str.length());
    if (!line_end) {
        line_end = str.end();
    }

    return {line_beg, line_end};
}

CString extract_line(CString& str) {
    const uint8* str_beg = str.ptr;
    const uint8* str_end = str.ptr + str.count;

    if (str_beg == str_end) {
        return {};
    }

    // If we start on a new line character or return character, skip these
    while (str_beg != str_end && (*str_beg == '\r' || *str_beg == '\n')) ++str_beg;

    const uint8* line_beg = str_beg;
    const uint8* line_end = (const uint8*)memchr(str_beg, '\n', str.length());
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
    const uint8* str_beg = str.data;
    const uint8* str_end = str.data + str.count;

    if (str_beg == str_end) {
        line = {};
        return false;
    }

    const uint8* line_beg = str_beg;
    const uint8* line_end = line_beg;

    // Find return or new line character
    while (line_end < str_end && *line_end != '\r' && *line_end != '\n') ++line_end;

    // Step over return or new line characters
    str_beg = MIN(line_end + 1, str_end);
    while (str_beg < str_end && (*str_beg == '\r' || *str_beg == '\n')) ++str_beg;

    // @NOTE: Do not modify line.count, its value contains the length of the buffer its pointing to
    auto count = MIN(line_end - line_beg, line.count - 1);
    line.data = (uint8*)memcpy(line.data, line_beg, count);
    line.data[count] = '\0';

    str.data = str_beg;
    str.count = str_end - str_beg;

    return true;
}
*/

ConversionResult<float32> to_float32(CString str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    float32 val = strtof(buf.cstr(), &end);
    return {val, end != buf.cstr()};
}

ConversionResult<float64> to_float64(CString str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    float64 val = strtod(buf.cstr(), &end);
    return {val, end != buf.cstr()};
}

ConversionResult<int32> to_int32(CString str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    int32 val = strtol(buf.cstr(), &end, 10);
    return {val, end != buf.cstr()};
}

ConversionResult<int64> to_int64(CString str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    int64 val = strtoll(buf.cstr(), &end, 10);
    return {val, end != buf.cstr()};
}

CString trim(CString str) {
    const uint8* beg = str.ptr;
    const uint8* end = str.ptr + str.count;

    while (beg < end && is_whitespace(*beg)) ++beg;
    while (end > beg && (is_whitespace(*(end - 1)) || *(end - 1) == '\0')) --end;

    return CString(beg, end - beg);
}

String trim(String str) {
    uint8* beg = str.ptr;
    uint8* end = str.ptr + str.count;

    while (beg < end && is_whitespace(*beg)) ++beg;
    while (end > beg && is_whitespace(*(end - 1))) --end;

    return String(beg, end - beg);
}

String allocate_and_read_textfile(CString filename) {
    StringBuffer<512> c_str_path = filename;
    FILE* file = fopen(c_str_path.cstr(), "rb");
    defer {
        if (file) fclose(file);
    };

    if (!file) return {};

// This is to handle big files. 64-bit versions
#ifdef _WIN32
#define FSEEK _fseeki64
#define FTELL _ftelli64
#else
#define FSEEK fseeko
#define FTELL ftello
#endif

    FSEEK(file, 0, SEEK_END);
    int64 file_size = FTELL(file);
    rewind(file);

    if (file_size <= 0) return {};

    uint8* ptr = (uint8*)MALLOC(file_size + 1);
    fread(ptr, 1, file_size, file);
    ptr[file_size] = '\0';

    return {ptr, file_size + 1};
}

CString get_directory(CString url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);

    const uint8* beg = url.begin();
    const uint8* end = url.end() - 1;

    while (end != beg && *end != '\\' && *end != '/') {
        end--;
    }

    return CString(beg, end - beg);
}

CString get_file(CString url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);

    const uint8* beg = url.end() - 1;
    const uint8* end = url.end();

    while (beg != url.begin() && *beg != '\\' && *beg != '/') {
        beg--;
    }
    if (*beg == '\\' || *beg == '/') beg++;

    return CString(beg, end - beg);
}

CString get_file_without_extension(CString url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);

    const uint8* beg = url.end() - 1;
    const uint8* end = url.end();

    while (beg != url.begin() && *beg != '\\' && *beg != '/') beg--;
    if (beg != url.begin()) beg++;

    while (end != beg && *end != '.') end--;

    return CString(beg, end - beg);
}

CString get_file_extension(CString url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);

    const uint8* beg = url.end() - 1;
    const uint8* end = url.end();

    while (beg != url.begin() && *beg != '.' && *beg != '\\' && *beg != '/') beg--;

    if (beg == url.begin() || *beg == '\\' || *beg == '/') {
        return CString();
    }

    beg++;  // skip '.'
    return CString(beg, end - beg);
}

inline static bool char_in_string(char c, CString str) {
    for (int64 i = 0; i < str.count; i++) {
        if (c == str[i]) return true;
    }
    return false;
}

StringBuffer<256> get_relative_path(CString from, CString to) {
    const uint8* c_from = from.beg();
    const uint8* c_to = to.beg();
    while (c_from != from.end() && c_to != to.end() && *c_from == *c_to) {
        c_from++;
        c_to++;
    }

    // If they have nothing in common. Return absolute path of to
    if (c_to == to.beg()) {
        return to;
    }

    int dir_count = 0;
    for (const uint8* c = c_from; c != from.end(); c++ /* <- LOL! */) {
        if (*c == '\\' || *c == '/') dir_count++;
    }

    StringBuffer<256> res;
    int offset = 0;
    for (int i = 0; i < dir_count; i++) {
        offset += snprintf(res.cstr() + offset, res.capacity(), "../");
    }

    StringBuffer<256> to_buf = CString(c_to, to.end());
    snprintf(res.cstr() + offset, res.capacity(), "%s", to_buf.beg());

    return res;
}

StringBuffer<256> get_absolute_path(CString absolute_reference, CString relative_file) {
    CString abs_dir = get_directory(absolute_reference);
    if (relative_file.count < 3) return {};

    StringBuffer<256> res;
    // If relative path is really an absolute path, just return that
    if (relative_file[0] == '/' || relative_file[1] == ':') {
        res = relative_file;
        return res;
    }

    int dir_count = 0;
    for (const uint8* c = relative_file.beg(); c < relative_file.end(); c += 3) {
        if (c[0] == '.' && (c + 1) != relative_file.end() && c[1] == '.' && (c + 2) != relative_file.end() && (c[2] == '/' || c[2] == '\\')) dir_count++;
    }

    const uint8* c = abs_dir.end() - 1;
    while (c > abs_dir.beg() && dir_count > 0) {
        if (*c == '/' || *c == '\\') {
            if (--dir_count == 0) break;
        }
        c--;
    }
    if (dir_count > 0 || c == abs_dir.beg()) return res;

    CString base_dir(abs_dir.beg(), c + 1);
    res = base_dir;
    StringBuffer<128> file = relative_file;
    snprintf(res.cstr() + base_dir.size(), res.capacity() - base_dir.size(), "/%s", file.beg());

    return res;
}

void convert_backslashes(String str) {
    for (uint8* c = str.beg(); c != str.end(); c++) {
        if (*c == '\\') *c = '/';
    }
}

bool is_digit(uint8 c) { return isdigit(c); }

bool is_alpha(uint8 c) { return isalpha(c); }

bool is_whitespace(uint8 c) { return isspace(c); }

bool contains_whitespace(CString str) {
    for (const uint8* c = str.beg(); c != str.end(); c++) {
        if (is_whitespace(*c)) return true;
    }
    return false;
}

bool balanced_parentheses(CString str) {
    int count = 0;
    const uint8* ptr = str.beg();
    while (ptr != str.end()) {
        if (*ptr == '(')
            count++;
        else if (*ptr == ')')
            count--;
        ptr++;
    }
    return count == 0;
}

CString extract_parentheses(CString str) {
    const uint8* beg = str.beg();

    while (beg != str.end() && *beg != '(') beg++;
    if (beg == str.end()) return {beg, str.end()};

    const uint8* end = beg + 1;
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

CString extract_parentheses_contents(CString str) {
    CString p = extract_parentheses(str);
    if (p.count < 2) return p;
    return {p.beg() + 1, p.end() - 1};
}

const uint8* find_character(CString str, uint8 c) { return (const uint8*)memchr(str.ptr, c, str.length()); }

bool contains_character(CString str, uint8 c) { return find_character(str, c) != str.end(); }

CString find_string(CString target, CString pattern) {
    if (target.count == 0 || pattern.count == 0) return {};

    char* ptr = Railgun_Trolldom((char*)target.ptr, (char*)pattern.ptr, (uint32)target.count, (uint32)pattern.count);
    if (ptr) {
        return {(uint8*)ptr, pattern.length()};
    }
    return {};
}

DynamicArray<String> tokenize(String str, uint8 delimiter) {
    DynamicArray<String> tokens;

    uint8* beg = str.beg();
    uint8* end = str.beg();

    while (end != str.end() && *end != '\0') {
        while (end != str.end() && *end != '\0' && *end != delimiter) end++;
        tokens.push_back(String(beg, end));
        beg = end;
        while (beg != str.end() && *end != '\0' && *beg == delimiter) beg++;
        end = beg;
    }

    return tokens;
}

DynamicArray<String> tokenize(String str, CString delimiter) {
    DynamicArray<String> tokens;

    uint8* beg = str.beg();
    uint8* end = str.beg();

    while (end != str.end() && *end != '\0') {
        while (end != str.end() && *end != '\0' && !char_in_string(*end, delimiter)) end++;
        tokens.push_back(String(beg, end));
        beg = end;
        while (beg != str.end() && *end != '\0' && char_in_string(*beg, delimiter)) beg++;
        end = beg;
    }

    return tokens;
}

DynamicArray<CString> ctokenize(CString str, uint8 delimiter) {
    DynamicArray<CString> tokens;

    const uint8* beg = str.beg();
    const uint8* end = str.beg();

    while (end != str.end() && *end != '\0') {
        while (end != str.end() && *end != '\0' && *end != delimiter) end++;
        tokens.push_back(CString(beg, end));
        beg = end;
        while (beg != str.end() && *beg != '\0' && *beg == delimiter) beg++;
        end = beg;
    }

    return tokens;
}

DynamicArray<CString> ctokenize(CString str, CString delimiter) {
    DynamicArray<CString> tokens;

    const uint8* beg = str.beg();
    const uint8* end = str.beg();

    while (end != str.end() && *end != '\0') {
        while (end != str.end() && *end != '\0' && !char_in_string(*end, delimiter)) end++;
        tokens.push_back(CString(beg, end));
        beg = end;
        while (beg != str.end() && *end != '\0' && char_in_string(*beg, delimiter)) beg++;
        end = beg;
    }

    return tokens;
}

constexpr uint8 delimiter = ':';
constexpr uint8 wildcard = '*';

// Range extraction functionality
bool is_range(CString arg) {
    for (const uint8* c = arg.beg(); c != arg.end(); c++) {
        if (is_digit(*c)) continue;
        if (*c == delimiter) return true;
        if (*c == wildcard) return true;
    }
    return false;
}

bool extract_range(IntRange* range, CString arg) {
    if (arg.count == 0) {
        *range = {-1, -1};
        return false;
    }

    if (arg.count == 1 && arg[0] == wildcard) {
        range->x = -1;
        range->y = -1;
        return true;
    }

    const uint8* mid = arg.beg();
    while (mid != arg.end() && *mid != delimiter) mid++;
    if (mid == arg.end()) return false;

    CString str_first(arg.beg(), mid);
    CString str_last(mid + 1, arg.end());

    if (str_first.count == 1 && str_first[0] == wildcard) {
        range->x = -1;
    } else {
        auto res = to_int32(str_first);
        if (!res) return false;
        range->x = res;
    }

    if (str_last.count == 1 && str_last[0] == wildcard) {
        range->y = -1;
    } else {
        auto res = to_int32(str_last);
        if (!res) return false;
        range->y = res;
    }

    return true;
}

bool extract_ranges(DynamicArray<IntRange>* ranges, Array<const CString> args) {
    ASSERT(ranges);

    for (auto arg : args) {
        if (is_range(arg)) {
            IntRange r;
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
