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

#if 0
// Blatantly stolen from libconfuse
// https://github.com/martinh/libconfuse/blob/master/src/fmemopen.c


/*
 * Copyright (c) 2017  Joachim Nilsson <troglobit@gmail.com>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>

FILE* fmemopen(void* buf, size_t len, const char* type) {
    int fd;
    FILE* fp;
    char tp[MAX_PATH - 13];
    char fn[MAX_PATH + 1];

    if (!GetTempPathA(sizeof(tp), tp)) return NULL;
    if (!GetTempFileNameA(tp, "confuse", 0, fn)) return NULL;

    fd = _open(fn, _O_CREAT | _O_RDWR | _O_SHORT_LIVED | _O_TEMPORARY | _O_BINARY, _S_IREAD | _S_IWRITE);
    if (fd == -1) return NULL;

    fp = _fdopen(fd, "w+");
    if (!fp) {
        _close(fd);
        return NULL;
    }

    fwrite(buf, len, 1, fp);
    rewind(fp);

    return fp;
}
#endif