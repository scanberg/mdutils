#pragma once

#include <core/string_types.h>
#include <core/string_utils.h>
#include <mol/aminoacid.h>

constexpr AminoAcid get_amino_acid_from_string(CStringView cstr) {

    // Skip leading numbers and crap
    while (cstr.count > 0 && !is_alpha(*cstr.beg())) {
        cstr.ptr++;
        cstr.count--;
    };

    if (cstr.count != 3) return AminoAcid::Unknown;
    const char seq[4] = {to_upper(cstr[0]), to_upper(cstr[1]), to_upper(cstr[2]), '\0'};
    CStringView sv = seq;

    for (unsigned int i = 0; i < aminoacid::num_amino_acids; i++) {
        CStringView sym = aminoacid::detail::symbols[i];
        if (sv == sym) return (AminoAcid)i;
    }
    return AminoAcid::Unknown;
}