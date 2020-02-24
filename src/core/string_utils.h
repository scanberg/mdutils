#pragma once

#include "string_types.h"

// Comparison of Strings
bool compare(CStringView str_a, CStringView str_b);
bool compare_ignore_case(CStringView str_a, CStringView str_b);

bool compare_n(CStringView str_a, CStringView str_b, i64 num_chars);
bool compare_n_ignore_case(CStringView str_a, CStringView str_b, i64 num_chars);

// Copy String
// Note: will zero terminate the dst String
void copy(StringView dst, CStringView src);

// Copy N first characters of src String
// Note: will zero terminate the dst String
void copy_n(StringView dst, CStringView src, i64 num_chars);

// Peeks into string and extracts the first line without modifying str.
CStringView peek_line(CStringView str);

// Reads the next line from a CString
// Note: Returns the extracted line, str gets modified.
// Note: NO guarantee that the line will be zero terminated. Be careful with printf!
CStringView extract_line(CStringView& str);

// Copies the next line from a CString
// Note: line holds the buffer which the line will be copied to, str gets modified
// Note: Guaranteed that the line will be zero terminated.
// String copy_line(String& line, CString& str);

StringView allocate_string(CStringView str);
StringView allocate_string(i32 length);
void free_string(StringView* str);

template <typename T>
struct ConversionResult {
    T value;
    bool success;
    constexpr operator T() const { return value; }
};

// Wrappers around strtof
ConversionResult<f32> to_float32(CStringView str);
ConversionResult<f64> to_float64(CStringView str);
inline ConversionResult<float> to_float(CStringView str) { return to_float32(str); }

// Wrappers around strtol
ConversionResult<i32> to_int32(CStringView str);
ConversionResult<i64> to_int64(CStringView str);
inline ConversionResult<int> to_int(CStringView str) { return to_int32(str); }

constexpr int char_to_digit(char c) { return c - '0'; }
constexpr char to_lower(char c) {
    if ('A' <= c && c <= 'Z') return c+32;
    return c;
}
constexpr char to_upper(char c) {
    if ('a' <= c && c <= 'z') return c-32;
    return c;
}
constexpr bool is_digit(char c) { return '0' <= c && c <= '9'; }
constexpr bool is_alpha(char c) { return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z'); }
constexpr bool is_whitespace(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }

template <typename Float = f32>
constexpr inline float str_to_float(CStringView str) {
	const char* c = str.beg();
    const char* end = str.end();

	Float val = 0;
	Float base = 1;
	Float sign = 1;

	while (c != end) {
		if ('0' <= *c && *c <= '9') {
			val = val * (Float)10 + char_to_digit(*c);
		}
		else if (*c == '-') {
			sign = -(Float)1;
		}
		else if (*c == '.') {
			c++;
			break;
		}
		c++;
	}
    while (c != end && *c != ' ') {
		if ('0' <= *c && *c <= '9') {
			base *= (Float)0.1;
			val += char_to_digit(*c) * base;
		}
		c++;
	}
    if (c != end && (*c == 'e' || *c == 'E')) {
		c++;
		// @TODO: Read exponent
	}

	return sign * val;
}

// Finds a character inside a string
constexpr const char* find_character(CStringView str, char character) {
    for (const char& c : str) {
        if (c == character) return &c;
    }
    // @FUTURE: Enable this if not evaluated at compile time
    // return (const char*)memchr(str.ptr, c, str.length());
    return nullptr;
}

// Finds a character inside a string
constexpr char* find_character(StringView str, char character) {
    for (char& c : str) {
        if (c == character) return &c;
    }
    // @FUTURE: Enable this if not evaluated at compile time
    // return (const char*)memchr(str.ptr, c, str.length());
    return nullptr;
}

constexpr bool contains_character(CStringView str, char c) { return find_character(str, c) != nullptr; }

constexpr bool contains_digits(CStringView str) {
    for (const char c : str) {
        if (is_digit(c)) return true;
    }
    return false;
}

constexpr bool contains_alpha(CStringView str) {
    for (const char c : str) {
        if (is_alpha(c)) return true;
    }
    return false;
}

constexpr bool contains_whitespace(CStringView str) {
    for (const char c : str) {
        if (is_whitespace(c)) return true;
    }
    return false;
}

constexpr bool balanced_parentheses(CStringView str) {
    i64 bal = 0;
    for (const char c : str) {
        if (c == '(') ++bal;
        else if (c == ')') --bal;
    }
    return bal == 0;
}

// Removes whitespace from begining and end of String
CStringView trim(CStringView str);
//String trim(String str);

// Reads text file and copies into allocated zero terminated String
StringView allocate_and_read_textfile(CStringView filename);

// Returns directory part from url, ex: func("C:/folder/file.ext") should return "C:/folder/"
CStringView get_directory(CStringView url);

// Returns file from url including extension, ex: func("C:/folder/file.ext") should return "file.ext"
CStringView get_file(CStringView url);

// Returns file part of url excluding extension, ex: func("C:/folder/file.ext") should return "file"
CStringView get_file_without_extension(CStringView url);

// Returns file extension part of url, ex: func("C:/folder/file.ext") should return "ext"
CStringView get_file_extension(CStringView url);

// Creates a relative path from the absolute path specified in the first file to the absolute path of the second file
StringBuffer<256> get_relative_path(CStringView absolute_from, CStringView absolute_to);

// Creates an absolute path using an absolute reference path specified in the first file and the relative path specified in the second
StringBuffer<256> get_absolute_path(CStringView absolute_reference, CStringView relative_file);

// Converts windows backslashes '\\' to forward slashes '/'
void convert_backslashes(StringView str);

CStringView extract_parentheses(CStringView str);
CStringView extract_parentheses_contents(CStringView str);

// Attempts to find a pattern inside a target string
// Returns empty CString if pattern is not found
// Otherwise it returns the CString pointing to the location of the pattern inside target and has the size of the pattern.
CStringView find_string(CStringView target, CStringView pattern);

// Tokenizes a string into shorter strings based on some delimiter
DynamicArray<CStringView> tokenize(CStringView str, char delimiter = ' ');
DynamicArray<CStringView> tokenize(CStringView str, CStringView delimiters);

// Positive range extraction functionality
// Examples of ranges:
// input -> output
// 1-4  -> { 1, 4}
// 1-*  -> { 1,-1}
// *-18 -> {-1,18}
// *-*  -> {-1,-1}

bool is_range(CStringView arg);
bool extract_range(Range<i32>* range, CStringView arg);
bool extract_ranges(DynamicArray<Range<i32>>* ranges, Array<const CStringView> args);

void print_string(CStringView str);