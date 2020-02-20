#include "file.h"

#if defined(_WIN32) && !defined(__CYGWIN__) && !defined(__GNUC__)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif
#include <windows.h>
#endif

FILE* fopen_utf8(const char* file, const char* mode) {

#if defined(_WIN32) && !defined(__CYGWIN__) && !defined(__GNUC__)
const int file_len = (int)strlen(file);
const int mode_len = (int)strlen(mode);

if (file_len == 0) return NULL;
if (mode_len == 0) return NULL;

wchar_t w_file[MAX_PATH];
wchar_t w_mode[MAX_PATH];

const int w_file_len = MultiByteToWideChar(CP_UTF8, 0, file, file_len, w_file, MAX_PATH);
if (w_file_len >= MAX_PATH) return NULL;
w_file[w_file_len] = L'\0';

const int w_mode_len = MultiByteToWideChar(CP_UTF8, 0, mode, mode_len, w_mode, MAX_PATH);
if (w_mode_len >= MAX_PATH) return NULL;
w_mode[w_mode_len] = L'\0';

return _wfopen(w_file, w_mode);

#else
    return fopen(file, mode);
#endif
}

FILE* fopen(CStringView file, CStringView mode) {
#if defined(_WIN32) && !defined(__CYGWIN__) && !defined(__GNUC__)
    if (file.length() == 0) return NULL;
    if (mode.length() == 0) return NULL;

    wchar_t w_file[MAX_PATH];
    wchar_t w_mode[MAX_PATH];

    const int w_file_len = MultiByteToWideChar(CP_UTF8, 0, file.cstr(), (int)file.length(), w_file, MAX_PATH);
    if (w_file_len >= MAX_PATH) return NULL;
    w_file[w_file_len] = L'\0';

    const int w_mode_len = MultiByteToWideChar(CP_UTF8, 0, mode.cstr(), (int)mode.length(), w_mode, MAX_PATH);
    if (w_mode_len >= MAX_PATH) return NULL;
    w_mode[w_mode_len] = L'\0';

    return _wfopen(w_file, w_mode);

#else
    StringBuffer<512> z_file = file;
    StringBuffer<512> z_mode = mode;
    return fopen(z_file.cstr(), z_mode.cstr());
#endif
}