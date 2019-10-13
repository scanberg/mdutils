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
using SeqIdx = int32;
using BondIdx = int32;

using AtomRange = Range<AtomIdx>;
using ResRange = Range<ResIdx>;
using ChainRange = Range<ChainIdx>;
using BondRange = Range<BondIdx>;

struct Bond {
    AtomIdx idx[2] = {0, 0};
};

struct BackboneSegment {
    AtomIdx ca_idx = -1;
    AtomIdx n_idx = -1;
    AtomIdx c_idx = -1;
    AtomIdx o_idx = -1;
};

struct BackboneAngle {
    float phi = 0;
    float psi = 0;
    operator vec2() const { return {phi, psi}; };
};

using BackboneSequence = ResRange;

struct HydrogenBondDonor {
    AtomIdx donor_idx = 0;
    AtomIdx hydro_idx = 0;
};

using HydrogenBondAcceptor = AtomIdx;

struct Residue {
    Label name{};
    ResIdx id = -1;
    // ChainIdx chain_idx = -1;
    // SeqIdx sequence_idx = -1;

    AtomRange atom_range;

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

struct Sequence {
    Label id = {};
    ResRange res_range = {};
    AtomRange atom_range = {};
};

using Chain = Sequence;

// Interface to access molecular data
struct MoleculeStructure {
    // SOA Layout for Atom data
    struct {
        // Aligned data
        int64 count = 0;
        struct {
            float* x = nullptr;
            float* y = nullptr;
            float* z = nullptr;
        } position;
        struct {
            float* x = nullptr;
            float* y = nullptr;
            float* z = nullptr;
        } velocity;
        float* radius = nullptr;
        float* mass = nullptr;

        // Non aligned data
        Element* element = nullptr;
        Label* label = nullptr;
        ResIdx* res_idx = nullptr;
        ChainIdx* chain_idx = nullptr;
        SeqIdx* seq_idx = nullptr;

    } atom;

    ArrayView<Bond> covalent_bonds{};
    ArrayView<Residue> residues{};
    ArrayView<Chain> chains{};
    ArrayView<Sequence> sequences{};

    struct {
        // Segments and angles should match in length and if not zero, they should match the number of residues
        ArrayView<BackboneSegment> segments{};
        ArrayView<BackboneAngle> angles{};
        ArrayView<BackboneSequence> sequences{};
    } backbone;

    struct {
        ArrayView<HydrogenBondDonor> donors{};
        ArrayView<HydrogenBondAcceptor> acceptors{};
    } hydrogen_bond;

    operator bool() const { return atom.count > 0; }
};

// General accessors
inline ArrayView<Residue> get_residues(MoleculeStructure& mol) { return mol.residues; }
inline ArrayView<const Residue> get_residues(const MoleculeStructure& mol) { return mol.residues; }
inline ArrayView<Chain> get_chains(MoleculeStructure& mol) { return mol.chains; }
inline ArrayView<const Chain> get_chains(const MoleculeStructure& mol) { return mol.chains; }
inline ArrayView<Bond> get_covalent_bonds(MoleculeStructure& mol) { return mol.covalent_bonds; }
inline ArrayView<const Bond> get_covalent_bonds(const MoleculeStructure& mol) { return mol.covalent_bonds; }

inline Float3Stream get_position_stream(MoleculeStructure& mol) { return {mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, mol.atom.count}; }
inline Float3Stream get_velocity_stream(MoleculeStructure& mol) { return {mol.atom.velocity.x, mol.atom.velocity.y, mol.atom.velocity.z, mol.atom.count}; }

// Single atom access
inline vec3 get_position_xyz(const MoleculeStructure& mol, AtomIdx idx) {
    ASSERT(0 <= idx && idx < mol.atom.count);
    return {mol.atom.position.x[idx], mol.atom.position.y[idx], mol.atom.position.z[idx]};
}

inline ArrayView<float> get_positions_x(MoleculeStructure& mol) { return ArrayView<float>(mol.atom.position.x, mol.atom.count); }
inline ArrayView<const float> get_positions_x(const MoleculeStructure& mol) { return ArrayView<const float>(mol.atom.position.x, mol.atom.count); }
inline ArrayView<float> get_positions_y(MoleculeStructure& mol) { return ArrayView<float>(mol.atom.position.y, mol.atom.count); }
inline ArrayView<const float> get_positions_y(const MoleculeStructure& mol) { return ArrayView<const float>(mol.atom.position.y, mol.atom.count); }
inline ArrayView<float> get_positions_z(MoleculeStructure& mol) { return ArrayView<float>(mol.atom.position.z, mol.atom.count); }
inline ArrayView<const float> get_positions_z(const MoleculeStructure& mol) { return ArrayView<const float>(mol.atom.position.z, mol.atom.count); }

inline ArrayView<float> get_velocities_x(MoleculeStructure& mol) { return ArrayView<float>(mol.atom.velocity.x, mol.atom.count); }
inline ArrayView<const float> get_velocities_x(const MoleculeStructure& mol) { return ArrayView<const float>(mol.atom.velocity.x, mol.atom.count); }
inline ArrayView<float> get_velocities_y(MoleculeStructure& mol) { return ArrayView<float>(mol.atom.velocity.y, mol.atom.count); }
inline ArrayView<const float> get_velocities_y(const MoleculeStructure& mol) { return ArrayView<const float>(mol.atom.velocity.y, mol.atom.count); }
inline ArrayView<float> get_velocities_z(MoleculeStructure& mol) { return ArrayView<float>(mol.atom.velocity.z, mol.atom.count); }
inline ArrayView<const float> get_velocities_z(const MoleculeStructure& mol) { return ArrayView<const float>(mol.atom.velocity.z, mol.atom.count); }

inline ArrayView<float> get_radii(MoleculeStructure& mol) { return ArrayView<float>(mol.atom.radius, mol.atom.count); }

inline ArrayView<Element> get_elements(MoleculeStructure& mol) { return ArrayView<Element>(mol.atom.element, mol.atom.count); }
inline ArrayView<const Element> get_elements(const MoleculeStructure& mol) { return ArrayView<const Element>(mol.atom.element, mol.atom.count); }
inline ArrayView<Label> get_labels(MoleculeStructure& mol) { return ArrayView<Label>(mol.atom.label, mol.atom.count); }
inline ArrayView<const Label> get_labels(const MoleculeStructure& mol) { return ArrayView<const Label>(mol.atom.label, mol.atom.count); }
inline ArrayView<ResIdx> get_residue_indices(MoleculeStructure& mol) { return ArrayView<ResIdx>(mol.atom.res_idx, mol.atom.count); }
inline ArrayView<const ResIdx> get_residue_indices(const MoleculeStructure& mol) { return ArrayView<const ResIdx>(mol.atom.res_idx, mol.atom.count); }

// Backbone accessors
inline ArrayView<BackboneSegment> get_backbone(MoleculeStructure& mol, BackboneSequence seq) {
    ASSERT(0 <= seq.beg && seq.end <= mol.backbone.segments.count);
    return mol.backbone.segments.subarray(seq);
}

inline ArrayView<const BackboneSegment> get_backbone(const MoleculeStructure& mol, BackboneSequence seq) {
    ASSERT(0 <= seq.beg && seq.end <= mol.backbone.segments.count);
    return mol.backbone.segments.subarray(seq);
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

inline Sequence get_sequence(MoleculeStructure& mol, SeqIdx idx) {
    ASSERT(0 <= idx && idx < mol.sequences.count);
    return mol.sequences[idx];
}

inline Sequence get_sequence(const MoleculeStructure& mol, SeqIdx idx) {
    ASSERT(0 <= idx && idx < mol.sequences.count);
    return mol.sequences[idx];
}

/*
inline Array<BackboneSegment> get_backbone(MoleculeStructure& mol, const Chain& chain) {
    ASSERT(0 <= chain.res_range.beg && chain.res_range.end <= mol.residues.count);
    if (mol.backbone.segments.count == 0) return {};
    return mol.backbone.segments.subarray(chain.res_range);
}

inline Array<const BackboneSegment> get_backbone(const MoleculeStructure& mol, const Chain& chain) {
    ASSERT(0 <= chain.res_range.beg && chain.res_range.end <= mol.residues.count);
    if (mol.backbone.segments.count == 0) return {};
    return mol.backbone.segments.subarray(chain.res_range);
}
*/

inline ArrayView<Residue> get_residues(MoleculeStructure& mol, const Chain& chain) { return mol.residues.subarray(chain.res_range); }
inline ArrayView<const Residue> get_residues(const MoleculeStructure& mol, const Chain& chain) { return mol.residues.subarray(chain.res_range); }

inline ArrayView<Element> get_elements(MoleculeStructure& mol, Chain& chain) { return get_elements(mol).subarray(chain.atom_range); }

inline ArrayView<const Element> get_elements(const MoleculeStructure& mol, const Chain& chain) { return get_elements(mol).subarray(chain.atom_range); }

inline ArrayView<Label> get_labels(MoleculeStructure& mol, Chain& chain) { return get_labels(mol).subarray(chain.atom_range); }

inline ArrayView<const Label> get_labels(const MoleculeStructure& mol, const Chain& chain) { return get_labels(mol).subarray(chain.atom_range); }

// Res func
inline ArrayView<Element> get_elements(MoleculeStructure& mol, const Residue& res) { return get_elements(mol).subarray(res.atom_range); }
inline ArrayView<const Element> get_elements(const MoleculeStructure& mol, const Residue& res) { return get_elements(mol).subarray(res.atom_range); }
inline ArrayView<Label> get_labels(MoleculeStructure& mol, const Residue& res) { return get_labels(mol).subarray(res.atom_range); }
inline ArrayView<const Label> get_labels(const MoleculeStructure& mol, const Residue& res) { return get_labels(mol).subarray(res.atom_range); }
inline ArrayView<Bond> get_bonds(MoleculeStructure& mol, const Residue& res) { return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg); }
inline ArrayView<const Bond> get_bonds(const MoleculeStructure& mol, const Residue& res) { return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg); }
inline ArrayView<Bond> get_internal_bonds(MoleculeStructure& mol, const Residue& res) { return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg); }
inline ArrayView<const Bond> get_internal_bonds(const MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.subarray(res.bond_idx.beg_internal, res.bond_idx.end_internal - res.bond_idx.beg_internal);
}

bool init_molecule_structure(MoleculeStructure* mol, int32 num_atoms, int32 num_bonds, int32 num_residues, int32 num_chains, int32 num_sequences, int32 num_backbone_segments = 0,
                             int32 num_backbone_sequences = 0, int32 num_hydrogen_bond_donors = 0, int32 num_hydrogen_bond_acceptors = 0);
void free_molecule_structure(MoleculeStructure* mol);
