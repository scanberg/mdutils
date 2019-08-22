#pragma once

#include "string_types.h"

// Comparison of Strings
bool compare(CString str_a, CString str_b);
bool compare_ignore_case(CString str_a, CString str_b);

bool compare_n(CString str_a, CString str_b, int64 num_chars);
bool compare_n_ignore_case(CString str_a, CString str_b, int64 num_chars);

// Copy String
// Note: will zero terminate the dst String
void copy(String dst, CString src);

// Copy N first characters of src String
// Note: will zero terminate the dst String
void copy_n(String dst, CString src, int64 num_chars);

// Peeks into string and extracts the first line without modifying str.
CString peek_line(CString str);

// Reads the next line from a CString
// Note: Returns the extracted line, str gets modified.
// Note: NO guarantee that the line will be zero terminated. Be careful with printf!
CString extract_line(CString& str);

// Copies the next line from a CString
// Note: line holds the buffer which the line will be copied to, str gets modified
// Note: Guaranteed that the line will be zero terminated.
// String copy_line(String& line, CString& str);

String allocate_string(CString str);
String allocate_string(int32 length);
void free_string(String* str);

template <typename T>
struct ConversionResult {
    T value;
    bool success;
    constexpr operator T() const { return value; }
};

// Wrappers around strtof
ConversionResult<float32> to_float32(CString str);
ConversionResult<float64> to_float64(CString str);
inline ConversionResult<float> to_float(CString str) { return to_float32(str); }

// Wrappers around strtol
ConversionResult<int32> to_int32(CString str);
ConversionResult<int64> to_int64(CString str);
inline ConversionResult<int> to_int(CString str) { return to_int32(str); }

inline int char_to_digit(char c) { return c - '0'; }

template <typename Float = float32>
constexpr inline float fast_str_to_float(CString str) {
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

// Removes whitespace from begining and end of String
CString trim(CString str);
//String trim(String str);

// Reads text file and copies into allocated zero terminated String
String allocate_and_read_textfile(CString filename);

// Returns directory part from url, ex: func("C:/folder/file.ext") should return "C:/folder/"
CString get_directory(CString url);

// Returns file from url including extension, ex: func("C:/folder/file.ext") should return "file.ext"
CString get_file(CString url);

// Returns file part of url excluding extension, ex: func("C:/folder/file.ext") should return "file"
CString get_file_without_extension(CString url);

// Returns file extension part of url, ex: func("C:/folder/file.ext") should return "ext"
CString get_file_extension(CString url);

// Creates a relative path from the absolute path specified in the first file to the absolute path of the second file
StringBuffer<256> get_relative_path(CString absolute_from, CString absolute_to);

// Creates an absolute path using an absolute reference path specified in the first file and the relative path specified in the second
StringBuffer<256> get_absolute_path(CString absolute_reference, CString relative_file);

// Converts windows backslashes '\\' to forward slashes '/'
void convert_backslashes(String str);

bool is_digit(char c);
bool is_alpha(char c);
bool is_whitespace(char c);
bool contains_whitespace(CString str);
bool balanced_parentheses(CString str);

CString extract_parentheses(CString str);
CString extract_parentheses_contents(CString str);

const char* find_character(CString str, char c);
bool contains_character(CString str, char c);

// Attempts to find a pattern inside a target string
// Returns empty CString if pattern is not found
// Otherwise it returns the CString pointing to the location of the pattern inside target and has the size of the pattern.
CString find_string(CString target, CString pattern);

// Tokenizes a string into shorter strings based on some delimiter
DynamicArray<CString> tokenize(CString str, char delimiter = ' ');
DynamicArray<CString> tokenize(CString str, CString delimiters);

// Positive range extraction functionality
// Examples of ranges:
// input -> output
// 1-4  -> { 1, 4}
// 1-*  -> { 1,-1}
// *-18 -> {-1,18}
// *-*  -> {-1,-1}

bool is_range(CString arg);
bool extract_range(Range<int32>* range, CString arg);
bool extract_ranges(DynamicArray<Range<int32>>* ranges, Array<const CString> args);

/*
// Temporary string object with managed memory
struct TmpString : CString {
    TmpString(CString str) {
        String tmp = allocate_string(str);
        this->ptr = tmp.ptr;
        this->count = tmp.count;
    }
    TmpString(const TmpString& other) = delete;
    TmpString(TmpString&& other) noexcept {
        this->ptr = other.ptr;
        this->count = other.count;
    }
    ~TmpString() { FREE((void*)this->ptr); }
};

// This is a hack to generate a zero terminated string from a CString object
// Returns an object with a temporary allocated string which is freed upon its destruction
inline TmpString make_tmp_str(CString str) { return TmpString(str); }
*/

inline void print_string(CString str) { printf("%.*s", (int)str.count, str.ptr); }

// TODO: Possibly implement a good templated print function in the style of printf as discussed here
// https://stackoverflow.com/questions/17671772/c11-variadic-printf-performance
