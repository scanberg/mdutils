#pragma once

#include "string_types.h"
#include <stdio.h>

// fopen with support utf-8 character encoding.

FILE* fopen_utf8(const char* filename, const char* mode);
FILE* fopen(CStringView filename, CStringView mode);