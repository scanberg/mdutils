#include <core/string_utils.h>
#include <core/common.h>
#include <stdio.h>
#include <ctype.h>

#ifdef WIN32
#pragma warning(disable : 4996)  // Unsafe strncpy
#endif

#define MIN(x, y) ((x < y) ? (x) : (y))
#define MAX(x, y) ((x > y) ? (x) : (y))

#include <stdint.h>
#include <stdlib.h>

#define ALPHABET_LEN 256
#define NOT_FOUND patlen

// delta1 table: delta1[c] contains the distance between the last
// character of pat and the rightmost occurrence of c in pat.
// If c does not occur in pat, then delta1[c] = patlen.
// If c is at string[i] and c != pat[patlen-1], we can
// safely shift i over by delta1[c], which is the minimum distance
// needed to shift pat forward to get string[i] lined up
// with some character in pat.
// this algorithm runs in alphabet_len+patlen time.
void make_delta1(int* delta1, uint8_t* pat, int32_t patlen) {
    int i;
    for (i = 0; i < ALPHABET_LEN; i++) {
        delta1[i] = NOT_FOUND;
    }
    for (i = 0; i < patlen - 1; i++) {
        delta1[pat[i]] = patlen - 1 - i;
    }
}

// true if the suffix of word starting from word[pos] is a prefix
// of word
int is_prefix(uint8_t* word, int wordlen, int pos) {
    int i;
    int suffixlen = wordlen - pos;
    // could also use the strncmp() library function here
    for (i = 0; i < suffixlen; i++) {
        if (word[i] != word[pos + i]) {
            return 0;
        }
    }
    return 1;
}

// length of the longest suffix of word ending on word[pos].
// suffix_length("dddbcabc", 8, 4) = 2
int suffix_length(uint8_t* word, int wordlen, int pos) {
    int i;
    // increment suffix length i to the first mismatch or beginning
    // of the word
    for (i = 0; (word[pos - i] == word[wordlen - 1 - i]) && (i < pos); i++)
        ;
    return i;
}

// delta2 table: given a mismatch at pat[pos], we want to align
// with the next possible full match could be based on what we
// know about pat[pos+1] to pat[patlen-1].
//
// In case 1:
// pat[pos+1] to pat[patlen-1] does not occur elsewhere in pat,
// the next plausible match starts at or after the mismatch.
// If, within the substring pat[pos+1 .. patlen-1], lies a prefix
// of pat, the next plausible match is here (if there are multiple
// prefixes in the substring, pick the longest). Otherwise, the
// next plausible match starts past the character aligned with
// pat[patlen-1].
//
// In case 2:
// pat[pos+1] to pat[patlen-1] does occur elsewhere in pat. The
// mismatch tells us that we are not looking at the end of a match.
// We may, however, be looking at the middle of a match.
//
// The first loop, which takes care of case 1, is analogous to
// the KMP table, adapted for a 'backwards' scan order with the
// additional restriction that the substrings it considers as
// potential prefixes are all suffixes. In the worst case scenario
// pat consists of the same letter repeated, so every suffix is
// a prefix. This loop alone is not sufficient, however:
// Suppose that pat is "ABYXCDBYX", and text is ".....ABYXCDEYX".
// We will match X, Y, and find B != E. There is no prefix of pat
// in the suffix "YX", so the first loop tells us to skip forward
// by 9 characters.
// Although superficially similar to the KMP table, the KMP table
// relies on information about the beginning of the partial match
// that the BM algorithm does not have.
//
// The second loop addresses case 2. Since suffix_length may not be
// unique, we want to take the minimum value, which will tell us
// how far away the closest potential match is.
void make_delta2(int* delta2, uint8_t* pat, int32_t patlen) {
    int p;
    int last_prefix_index = patlen - 1;

    // first loop
    for (p = patlen - 1; p >= 0; p--) {
        if (is_prefix(pat, patlen, p + 1)) {
            last_prefix_index = p + 1;
        }
        delta2[p] = last_prefix_index + (patlen - 1 - p);
    }

    // second loop
    for (p = 0; p < patlen - 1; p++) {
        int slen = suffix_length(pat, patlen, p);
        if (pat[p - slen] != pat[patlen - 1 - slen]) {
            delta2[patlen - 1 - slen] = patlen - 1 - p + slen;
        }
    }
}

// This is taken from https://en.wikipedia.org/wiki/Boyer%E2%80%93Moore_string-search_algorithm
uint8_t* boyer_moore(uint8_t* string, uint32_t stringlen, uint8_t* pat, uint32_t patlen) {
    int i;
    int delta1[ALPHABET_LEN];
    int* delta2 = (int*)malloc(patlen * sizeof(int));
    make_delta1(delta1, pat, patlen);
    make_delta2(delta2, pat, patlen);

    // The empty pattern must be considered specially
    if (patlen == 0) {
        free(delta2);
        return string;
    }

    i = patlen - 1;
    while (i < stringlen) {
        int j = patlen - 1;
        while (j >= 0 && (string[i] == pat[j])) {
            --i;
            --j;
        }
        if (j < 0) {
            free(delta2);
            return (string + i + 1);
        }

        i += MAX(delta1[string[i]], delta2[j]);
    }
    free(delta2);
    return NULL;
}

static inline bool internal_compare(const uint8* str_a, const uint8* str_b, int64 len, bool ignore_case) {
    if (ignore_case) {
        for (int64 i = 0; i < len; i++) {
            if (tolower(str_a[i]) != tolower(str_b[i])) return false;
        }
    } else {
        for (int64 i = 0; i < len; i++) {
            if (str_a[i] != str_b[i]) return false;
        }
    }
    return true;
}

bool compare(CString str_a, CString str_b, bool ignore_case) {
    // int64 len = MIN(str_a.count, str_b.count);
    if (str_a.count != str_b.count) return false;
    if (str_a.count == 0) return false;
    return internal_compare(str_a.data, str_b.data, str_a.count, ignore_case);
}

bool compare_n(CString str_a, CString str_b, int64 num_chars, bool ignore_case) {
    int64 len = MIN(str_a.count, MIN(str_b.count, num_chars));
    if (len < num_chars) return false;
    return internal_compare(str_a.data, str_b.data, len, ignore_case);
}

String find_needle(CString needle, String haystack) {
    uint8* res = boyer_moore(haystack.data, haystack.count, (uint8*)needle.data, needle.count);
    if (res) return {res, needle.length()};
    return {};
}

CString find_needle(CString needle, CString haystack) {
	uint8* res = boyer_moore((uint8*)haystack.data, haystack.count, (uint8*)needle.data, needle.count);
    if (res) return {res, needle.length()};
    return {};
}

void copy(String dst, CString src) {
    ASSERT(dst.data != 0);
    ASSERT(src.data != 0);
    auto len = MIN(dst.count, src.count);
    memcpy(dst.data, src.data, len);
    dst.data[src.count] = '\0';
}

void copy_n(String dst, CString src, int64 num_chars) {
    ASSERT(dst.data != 0);
    ASSERT(src.data != 0);
    auto len = MIN(dst.count, src.count);
    len = MIN(len, num_chars);
    memcpy(dst.data, src.data, len);
    dst.data[num_chars] = '\0';
}

String allocate_string(CString str) {
    if (str.count == 0) return {};
    uint8* data = (uint8*)MALLOC(str.count + 1);
    strncpy((char*)data, str.cstr(), str.count + 1);
    return {data, str.count};
}

String allocate_string(int32 length) {
    if (length == 0) return {};
    uint8* data = (uint8*)MALLOC(length);
    return {data, length};
}

void free_string(String* str) {
    if (str->data) {
        FREE(str->data);
        str->data = nullptr;
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
    const uint8* str_beg = str.data;
    const uint8* str_end = str.data + str.count;

    if (str_beg == str_end) {
        return {};
    }

    const uint8* line_beg = str_beg;
    const uint8* line_end = (const uint8*)memchr(str_beg, '\n', str.length());
    if (!line_end) {
        line_end = str_end;
        str_beg = str_end;
    }
	else {
        // Step over return and new line characters
        str_beg = MIN(line_end + 1, str_end);
        while (str_beg != str_end && (*str_beg == '\r' || *str_beg == '\n')) ++str_beg;
	}

    str.data = str_beg;
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
    float32 val = strtof(buf, &end);
    return {val, end != buf.cstr()};
}

ConversionResult<float64> to_float64(CString str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    float64 val = strtod(buf, &end);
    return {val, end != buf.cstr()};
}

ConversionResult<int32> to_int32(CString str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    int32 val = strtol(buf, &end, 10);
    return {val, end != buf.cstr()};
}

ConversionResult<int64> to_int64(CString str) {
    // Make sure that the string passed into atof is zero-terminated
    StringBuffer<32> buf = str;
    char* end = nullptr;
    int64 val = strtoll(buf, &end, 10);
    return {val, end != buf.cstr()};
}

CString trim(CString str) {
    const uint8* beg = str.data;
    const uint8* end = str.data + str.count;

    while (beg < end && is_whitespace(*beg)) ++beg;
    while (end > beg && (is_whitespace(*(end - 1)) || *(end - 1) == '\0')) --end;

    return CString(beg, end - beg);
}

String trim(String str) {
    uint8* beg = str.data;
    uint8* end = str.data + str.count;

    while (beg < end && is_whitespace(*beg)) ++beg;
    while (end > beg && is_whitespace(*(end - 1))) --end;

    return String(beg, end - beg);
}

String allocate_and_read_textfile(CString filename) {
    StringBuffer<512> c_str_path = filename;
    FILE* file = fopen(c_str_path.cstr(), "rb");
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

    uint8* data = (uint8*)MALLOC(file_size + 1);
    fread(data, 1, file_size, file);
    data[file_size] = '\0';

    return {data, file_size + 1};
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

    return CString(beg + 1, end - beg);
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
        offset += snprintf(res.cstr() + offset, res.MAX_LENGTH, "../");
    }

    StringBuffer<256> to_buf = CString(c_to, to.end());
    snprintf(res.cstr() + offset, res.MAX_LENGTH, "%s", to_buf.beg());

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
    StringBuffer<128> file = get_file(relative_file);
    snprintf(res.cstr() + base_dir.count, res.MAX_LENGTH - base_dir.count, "/%s", file.beg());

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

const uint8* find_character(CString str, uint8 c) {
    const uint8* ptr = str.beg();
    while (ptr < str.end() && *ptr != c) ptr++;
    return ptr;
}

bool contains_character(CString str, uint8 c) { return find_character(str, c) != str.end(); }

CString find_first_match(CString str, CString match) {
    if (str.count == 0 || match.count == 0) return {};

    const uint8* ptr = str.beg();

    while (ptr != str.end()) {
        if (*ptr == *match.beg()) {
            CString candidate(ptr, MIN(str.end() - ptr, match.count));
            if (compare(candidate, match)) return candidate;
        }
        ptr++;
    }
    return {};
}

bool contains_string(CString big_str, CString str) { return (bool)find_first_match(big_str, str); }

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
