#pragma once

#include <core/string_types.h>
#include <core/string_utils.h>
#include <mol/element.h>

constexpr Element get_element_from_string(CStringView cstr, bool ignore_case = false) {
    if (cstr.count == 0) return Element::Unknown;

    // Prune
    const auto* c = cstr.beg();
    while (c != cstr.end() && !is_alpha(*c)) c++;
    if (c == cstr.end()) {
        return Element::Unknown;
    }
    cstr.ptr = c;
    while (c != cstr.end() && is_alpha(*c)) c++;
    cstr.count = c - cstr.beg();

    // Two or more (Try to match first two)
    if (ignore_case) {
        if (cstr.length() > 1) {
            for (int32 i = 0; i < element::num_elements; i++) {
                CStringView elem = element::symbol((Element)i);
                if (elem.size() == 2 && cstr[0] == elem[0] && to_lower(cstr[1]) == elem[1]) return (Element)i;
            }
        }
    } else {
        if (cstr.length() > 1) {
            for (int32 i = 0; i < element::num_elements; i++) {
                CStringView elem = element::symbol((Element)i);
                if (elem.size() == 2 && cstr[0] == elem[0] && cstr[1] == elem[1]) return (Element)i;
            }
        }
    }

    // Try to match against first character
    for (int32 i = 0; i < element::num_elements; i++) {
        CStringView elem = element::symbol((Element)i);
        if (elem.size() == 1 && cstr[0] == elem[0]) return (Element)i;
    }

    return Element::Unknown;
}