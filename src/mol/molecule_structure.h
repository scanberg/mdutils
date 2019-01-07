#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/string_utils.h>
#include <mol/element.h>

using Label = StringBuffer<8>;
using AtomIdx = int32;
using ResIdx = int32;
using ChainIdx = int32;
using BondIdx = int32;

struct Bond {
    AtomIdx idx[2] = {0, 0};
};

struct BackboneSegment {
    AtomIdx ca_idx = -1;
    AtomIdx n_idx = -1;
    AtomIdx c_idx = -1;
    AtomIdx o_idx = -1;
};

struct BackboneSequence {
    ResIdx beg = 0;
    ResIdx end = 0;
};

struct HydrogenBondDonor {
    AtomIdx donor_idx = 0;
    AtomIdx hydro_idx = 0;
};

typedef AtomIdx HydrogenBondAcceptor;

struct Residue {
    Label name{};
    ResIdx id = -1;
    ChainIdx chain_idx = -1;

    struct {
        AtomIdx beg = 0;
        AtomIdx end = 0;
    } atom_idx;

    struct {
        // Covalent bonds for a residue
        // [beg, end[						the range for all bonds connected to this residue.
        // [beg, beg_internal[				the range for external shared bonds to previous residue in the chain.
        // [beg_internal, end_internal[		the range for internal bonds.
        // [end_internal, end[				the range for external shared bonds to next residue in the chain.
        BondIdx beg = 0;
        BondIdx beg_internal = 0;
        BondIdx end_internal = 0;
        BondIdx end = 0;
    } bond_idx;
};

struct Chain {
    Label id{};
    struct {
        ResIdx beg = 0;
        ResIdx end = 0;
    } res_idx;
};

// Interface to access molecular data
struct MoleculeStructure {
    struct {
        int64 count = 0;
        vec3* positions = nullptr;
        Element* elements = nullptr;
        Label* labels = nullptr;
        ResIdx* residue_indices = nullptr;
    } atom;

    Array<Bond> covalent_bonds{};
    Array<Residue> residues{};
    Array<Chain> chains{};

    // If this is not zero in length it should have the same length as residues
    Array<BackboneSegment> backbone_segments{};
    Array<BackboneSequence> backbone_sequences{};

    struct {
        Array<HydrogenBondDonor> donors{};
        Array<HydrogenBondAcceptor> acceptors{};
    } hydrogen_bond;

    operator bool() const { return atom.count > 0; }
};

// Atom data accessors
inline Array<vec3> get_positions(MoleculeStructure& mol) { return Array<vec3>(mol.atom.positions, mol.atom.count); }
inline Array<const vec3> get_positions(const MoleculeStructure& mol) { return Array<const vec3>(mol.atom.positions, mol.atom.count); }
inline Array<Element> get_elements(MoleculeStructure& mol) { return Array<Element>(mol.atom.elements, mol.atom.count); }
inline Array<const Element> get_elements(const MoleculeStructure& mol) { return Array<const Element>(mol.atom.elements, mol.atom.count); }
inline Array<Label> get_labels(MoleculeStructure& mol) { return Array<Label>(mol.atom.labels, mol.atom.count); }
inline Array<const Label> get_labels(const MoleculeStructure& mol) { return Array<const Label>(mol.atom.labels, mol.atom.count); }
inline Array<ResIdx> get_residue_indices(MoleculeStructure& mol) { return Array<ResIdx>(mol.atom.residue_indices, mol.atom.count); }
inline Array<const ResIdx> get_residue_indices(const MoleculeStructure& mol) { return Array<const ResIdx>(mol.atom.residue_indices, mol.atom.count); }

// Backbone accessors
inline Array<BackboneSegment> get_backbone(MoleculeStructure& mol, BackboneSequence seq) {
    ASSERT(0 <= seq.beg && seq.end <= mol.backbone_segments.count);
    return mol.backbone_segments.sub_array(seq.beg, seq.end - seq.beg);
}

inline Array<const BackboneSegment> get_backbone(const MoleculeStructure& mol, BackboneSequence seq) {
    ASSERT(0 <= seq.beg && seq.end <= mol.backbone_segments.count);
    return mol.backbone_segments.sub_array(seq.beg, seq.end - seq.beg);
}

// Chain accessors
inline Chain get_chain(MoleculeStructure& mol, ChainIdx idx) {
    ASSERT(0 <= idx && idx < mol.chains.count);
    return mol.chains[idx];
}

inline Chain get_chain(const MoleculeStructure& mol, ChainIdx idx) {
    ASSERT(0 <= idx && idx < mol.chains.count);
    return mol.chains[idx];
}

inline Array<BackboneSegment> get_backbone(MoleculeStructure& mol, Chain chain) {
    ASSERT(0 <= chain.res_idx.beg && chain.res_idx.end <= mol.residues.count);
    if (mol.backbone_segments.count == 0) return {};
    return mol.backbone_segments.sub_array(chain.res_idx.beg, chain.res_idx.end - chain.res_idx.beg);
}

inline Array<const BackboneSegment> get_backbone(const MoleculeStructure& mol, Chain chain) {
    ASSERT(0 <= chain.res_idx.beg && chain.res_idx.end <= mol.residues.count);
    if (mol.backbone_segments.count == 0) return {};
    return mol.backbone_segments.sub_array(chain.res_idx.beg, chain.res_idx.end - chain.res_idx.beg);
}

inline Array<Residue> get_residues(MoleculeStructure& mol, Chain chain) {
    return mol.residues.sub_array(chain.res_idx.beg, chain.res_idx.end - chain.res_idx.beg);
}

inline Array<const Residue> get_residues(const MoleculeStructure& mol, Chain chain) {
    return mol.residues.sub_array(chain.res_idx.beg, chain.res_idx.end - chain.res_idx.beg);
}

inline AtomIdx get_atom_beg_idx(MoleculeStructure& mol, Chain chain) {
    ASSERT(0 <= chain.res_idx.beg && chain.res_idx.beg < mol.residues.count);
    return mol.residues[chain.res_idx.beg].atom_idx.beg;
}

inline AtomIdx get_atom_beg_idx(const MoleculeStructure& mol, Chain chain) {
    ASSERT(0 <= chain.res_idx.beg && chain.res_idx.beg < mol.residues.count);
    return mol.residues[chain.res_idx.beg].atom_idx.beg;
}

inline AtomIdx get_atom_end_idx(MoleculeStructure& mol, Chain chain) {
    ASSERT(0 < chain.res_idx.end && chain.res_idx.end <= mol.residues.count);
    return mol.residues[chain.res_idx.end - 1].atom_idx.end;
}

inline AtomIdx get_atom_end_idx(const MoleculeStructure& mol, Chain chain) {
    ASSERT(0 < chain.res_idx.end && chain.res_idx.end <= mol.residues.count);
    return mol.residues[chain.res_idx.end - 1].atom_idx.end;
}

inline Array<vec3> get_positions(MoleculeStructure& mol, Chain chain) {
    const auto beg = get_atom_beg_idx(mol, chain);
    const auto end = get_atom_end_idx(mol, chain);
    return get_positions(mol).sub_array(beg, end - beg);
}

inline Array<const vec3> get_positions(const MoleculeStructure& mol, Chain chain) {
    const auto beg = get_atom_beg_idx(mol, chain);
    const auto end = get_atom_end_idx(mol, chain);
    return get_positions(mol).sub_array(beg, end - beg);
}

inline Array<Element> get_elements(MoleculeStructure& mol, Chain chain) {
    const auto beg = get_atom_beg_idx(mol, chain);
    const auto end = get_atom_end_idx(mol, chain);
    return get_elements(mol).sub_array(beg, end - beg);
}

inline Array<const Element> get_elements(const MoleculeStructure& mol, Chain chain) {
    const auto beg = get_atom_beg_idx(mol, chain);
    const auto end = get_atom_end_idx(mol, chain);
    return get_elements(mol).sub_array(beg, end - beg);
}

inline Array<Label> get_labels(MoleculeStructure& mol, Chain chain) {
    const auto beg = get_atom_beg_idx(mol, chain);
    const auto end = get_atom_end_idx(mol, chain);
    return get_labels(mol).sub_array(beg, end - beg);
}

inline Array<const Label> get_labels(const MoleculeStructure& mol, Chain chain) {
    const auto beg = get_atom_beg_idx(mol, chain);
    const auto end = get_atom_end_idx(mol, chain);
    return get_labels(mol).sub_array(beg, end - beg);
}

// Res func
inline Array<vec3> get_positions(MoleculeStructure& mol, const Residue& res) {
    return get_positions(mol).sub_array(res.atom_idx.beg, res.atom_idx.end - res.atom_idx.beg);
}

inline Array<const vec3> get_positions(const MoleculeStructure& mol, const Residue& res) {
    return get_positions(mol).sub_array(res.atom_idx.beg, res.atom_idx.end - res.atom_idx.beg);
}

inline Array<Element> get_elements(MoleculeStructure& mol, const Residue& res) {
    return get_elements(mol).sub_array(res.atom_idx.beg, res.atom_idx.end - res.atom_idx.beg);
}

inline Array<const Element> get_elements(const MoleculeStructure& mol, const Residue& res) {
    return get_elements(mol).sub_array(res.atom_idx.beg, res.atom_idx.end - res.atom_idx.beg);
}

inline Array<Label> get_labels(MoleculeStructure& mol, const Residue& res) {
    return get_labels(mol).sub_array(res.atom_idx.beg, res.atom_idx.end - res.atom_idx.beg);
}

inline Array<const Label> get_labels(const MoleculeStructure& mol, const Residue& res) {
    return get_labels(mol).sub_array(res.atom_idx.beg, res.atom_idx.end - res.atom_idx.beg);
}

inline Array<Bond> get_bonds(MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.sub_array(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg);
}

inline Array<const Bond> get_bonds(const MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.sub_array(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg);
}

inline Array<Bond> get_internal_bonds(MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.sub_array(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg);
}

inline Array<const Bond> get_internal_bonds(const MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.sub_array(res.bond_idx.beg_internal, res.bond_idx.end_internal - res.bond_idx.beg_internal);
}

bool init_molecule_structure(MoleculeStructure* mol, int32 num_atoms, int32 num_bonds, int32 num_residues, int32 num_chains,
                             int32 num_backbone_segments = 0, int32 num_backbone_sequences = 0, int32 num_hydrogen_bond_donors = 0,
                             int32 num_hydrogen_bond_acceptors = 0);
void free_molecule_structure(MoleculeStructure* mol);
