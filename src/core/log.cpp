#include <core/log.h>
#include <stdarg.h>
#include <stdio.h>
#include <time.h>

namespace logging {

struct LoggingBackend {
    LoggingFunc func = nullptr;
    void* usr_data = nullptr;
};

static DynamicArray<LoggingBackend> entries;
static const int buf_size = KILOBYTES(64); // @NOTE (Robin): 64K Ought to be enough for everyone, right?!
static char buf[buf_size];

void initialize() {}
void shutdown() {}

void register_backend(LoggingFunc func, void* usr_data) { entries.push_back({func, usr_data}); }

void record(Severity severity, const char* format, ...) {
    time_t now = time(0);
    struct tm tstruct;
    tstruct = *localtime(&now);

    va_list ap;
    va_start(ap, format);
    int count = 0;
    count += (int)strftime(buf, buf_size, "%T: ", &tstruct);
    count += vsnprintf(buf + count, buf_size, format, ap);
    va_end(ap);
    CStringView str(buf, count);

    for (const auto& e : entries) {
        e.func(str, severity, e.usr_data);
    }
}

}  // namespace logging