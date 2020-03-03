#pragma once

#define PLATFORM_WINDOWS 0
#define PLATFORM_LINUX   0
#define PLATFORM_OSX     0

#if defined _WIN32 || defined _WIN64
#ifndef NOMINMAX
#define NOMINMAX
#endif  // NOMINMAX

#undef  PLATFORM_WINDOWS
#define PLATFORM_WINDOWS 1

#elif __APPLE__
#include "TargetConditionals.h"
#ifdef TARGET_OS_MAC
#undef  PLATFORM_OSX
#define PLATFORM_OSX 1
#endif

#elif defined __linux__
#undef  PLATFORM_LINUX
#define PLATFORM_LINUX 1
#else

#error "Platform Unsuported"

#endif