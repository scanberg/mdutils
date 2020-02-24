#pragma once

#include "string_types.h"
#include <stdio.h>

// fopen with support utf-8 character encoding.

FILE* fopen_utf8(const char* filename, const char* mode);
FILE* fopen(CStringView filename, CStringView mode);

int64_t ftelli64(FILE* file);
int fseeki64(FILE* file, int64_t offset, int origin);