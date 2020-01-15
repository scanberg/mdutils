#pragma once

// clang-format off
enum class AminoAcid : unsigned char {
	Unknown = 0,
	Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly,
	His, Ile, Leu, Lys, Met, Phe, Pro, Ser,
	Thr, Trp, Tyr, Val, Sec, Pyl, Asx, Glx,
	Xle
};
// clang-format on

namespace aminoacid {
constexpr int num_amino_acids = 27;

namespace detail {
constexpr const char* symbols[num_amino_acids] = {"XAA", "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU",
                                                         "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                                                         "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASC", "GLX", "XLE"};
constexpr const char* names[num_amino_acids] = {"Unknown",
                                                "Alanine",
                                                "Arginine",
                                                "Asparagine",
                                                "Aspartic acid",
                                                "Cysteine",
                                                "Cysteine XX",
                                                "Glutamic acid",
                                                "Glutamine",
                                                "Glycine",
                                                "Histidine",
                                                "Isoleucine",
                                                "Luecine",
                                                "Lysine",
                                                "Methionine",
                                                "Phenylalanine",
                                                "Proline",
                                                "Serine",
                                                "Threonine",
                                                "Tryptophan",
                                                "Tyrosine",
                                                "Valine",
                                                "Selencysteine",
                                                "Pyrrolysine",
                                                "Asparagine or Aspartic acid",
                                                "Glutamine or glutamic acid",
                                                "Leucine or Isloleucine"};
}  // namespace detail

constexpr const char* name(AminoAcid amino) { return detail::names[(int)amino]; }
constexpr const char* symbol(AminoAcid amino) { return detail::symbols[(int)amino]; }
// @TODO: Add hydrophobicity and other animo-acid parameters
}  // namespace aminoacid