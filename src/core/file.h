#pragma once

#include "string_types.h"
#include <stdio.h>

// Wrapper for fopen to support utf-8 character encoding.

FILE* fopen_utf8(const char* filename, const char* mode);
FILE* fopen(CStringView filename, CStringView mode);

FILE* fmemopen(void* buf, size_t len, const char* mode);