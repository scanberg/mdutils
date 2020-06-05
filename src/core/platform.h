#pragma once

#define PLATFORM_WINDOWS 0
#define PLATFORM_LINUX   0
#define PLATFORM_OSX     0

#define COMPILER_GCC     0
#define COMPILER_CLANG   0
#define COMPILER_MSVC    0
#define COMPILER_UNKNOWN 0

#if defined(_WIN32) || defined(_WIN64)
    #undef  PLATFORM_WINDOWS
    #define PLATFORM_WINDOWS 1
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
#elif __APPLE__
    #include "TargetConditionals.h"
    #ifdef  TARGET_OS_MAC
        #undef  PLATFORM_OSX
        #define PLATFORM_OSX 1
    #endif
#elif defined __linux__
    #undef  PLATFORM_LINUX
    #define PLATFORM_LINUX 1
#else
    #error "Platform Unsuported"
#endif

#if defined(_MSC_VER)
#undef  COMPILER_MSVC
#define COMPILER_MSVC 1
#elif defined(__GNUC__) && !defined(__clang__)
    #undef  COMPILER_GCC
    #define COMPILER_GCC 1
#elif defined(__clang__)
    #undef  COMPILER_CLANG
    #define COMPILER_CLANG 1
#else
    #undef  COMPILER_UNKNOWN
    #define COMPILER_UNKNOWN 1
#endif