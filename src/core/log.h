#pragma once

#include "string_utils.h"

#define LOG_NOTE(...) logging::record(logging::Severity::Note, __VA_ARGS__);
#define LOG_WARNING(...) logging::record(logging::Severity::Warning, __VA_ARGS__);
#define LOG_ERROR(...) logging::record(logging::Severity::Error, __VA_ARGS__);
#define LOG_FATAL(...) logging::record(logging::Severity::Fatal, __VA_ARGS__);

namespace logging {
enum class Severity { Note, Warning, Error, Fatal };

typedef void (*LoggingFunc)(CStringView str, Severity severity, void* usr_data);

void initialize();
void shutdown();

void register_backend(LoggingFunc func, void* usr_data = nullptr);

void record(Severity severity, const char* format, ...);
}  // namespace logging
