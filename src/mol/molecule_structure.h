#pragma once

#include <core/types.h>
#include <core/vector_types.h>
#include <core/string_utils.h>
#include <mol/element.h>

using Label = StringBuffer<8>;
using AtomIdx = i32;
using ResIdx = i32;
using ChainIdx = i32;
using SegIdx = i32;
using BondIdx = i32;

using AtomRange = Range<AtomIdx>;
using ResRange = Range<ResIdx>;
using SegRange = Range<SegIdx>;
using BondRange = Range<BondIdx>;

using BackboneSequence = ResRange;
using HydrogenBondAcceptor = AtomIdx;

enum class SecondaryStructure : u8 { Undefined, Coil = Undefined, Sheet, Helix };

struct Bond {
    AtomIdx idx[2] = {0, 0};
};

struct BackboneSegment {
    AtomIdx ca_idx = 0;
    AtomIdx n_idx = 0;
    AtomIdx c_idx = 0;
    AtomIdx o_idx = 0;
};

struct BackboneAngle {
    float phi = 0;
    float psi = 0;
    operator vec2() const { return {phi, psi}; };
};

struct HydrogenBondDonor {
    AtomIdx donor_idx = 0;
    AtomIdx hydro_idx = 0;
};

/*
struct Residue {
    Label name{};
    ResIdx id = -1;
    // ChainIdx chain_idx = -1;
    // SeqIdx sequence_idx = -1;

    AtomRange atom_range{};

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
    ResRange res_range = {};
    AtomRange atom_range = {};
};

struct Chain : Sequence {
    Label id{};
};
*/

// Interface to access molecular data
struct MoleculeStructure {
    // SOA Layout for atom data
    struct {
        // Aligned data
        i64 count = 0;
        struct {
            float* x = nullptr;
            float* y = nullptr;
            float* z = nullptr;
        } position;
        float* radius = nullptr;
        float* mass = nullptr;

        // Non aligned data
        Element* element = nullptr;
        Label* label = nullptr;
        ResIdx* res_idx = nullptr;
        ChainIdx* chain_idx = nullptr;
        //SeqIdx* seq_idx = nullptr;
    } atom;

    struct {
        i64 count = 0;
        Bond* bond = nullptr;
    } covalent_bond;

    /*
    Array<Bond> covalent_bonds{};
    Array<Residue> residues{};
    Array<Chain> chains{};
    Array<Sequence> sequences{};
    */

    // SOA layout for residue data
    struct {
        i64 count = 0;
        ResIdx* id = nullptr;
        Label* name = nullptr;
        AtomRange* atom_range = nullptr;
        struct {
            // External bonds are implicitly given by the range in complete but not in internal
            // Bonds to previous residue in chain is given by: [complete.beg, internal.beg]
            // Bonds to next residue in chain is given by: [internal.end, complete.end]
            BondRange* complete = nullptr;
            BondRange* intra = nullptr;
        } bond;
    } residue;

    // Chains represent connected residues (through covalent bonds) and are usually defined as a single character starting from 'A'
    struct {
        i64 count = 0;
        Label* id = nullptr;
        ResRange* residue_range = nullptr;
    } chain;

    // Backbone represents
    struct {
        struct {
            i64 count = 0;
            BackboneSegment* segment = nullptr;
            BackboneAngle* angle = nullptr;
            SecondaryStructure* secondary_structure = nullptr;
        } segment;

        struct {
            i64 count = 0;
            SegRange* segment_range = nullptr;
        } sequence;

        /*
        // Segments and angles should match in length and if not zero, they should match the number of residues
        Array<BackboneSegment> segments{};
        Array<SecondaryStructure> secondary_structures{};
        Array<BackboneAngle> angles{};
        Array<BackboneSequence> sequences{};
        */
    } backbone;

    struct {
        struct {
            i64 count = 0;
            HydrogenBondDonor* data = nullptr;
        } donor;
        struct {
            i64 count = 0;
            HydrogenBondAcceptor* data = nullptr;
        } acceptor;
    } hydrogen_bond;

    operator bool() const { return atom.count > 0; }
};

inline vec3 get_atom_position(const MoleculeStructure& mol, AtomIdx i) {
    ASSERT(0 <= i && i <= mol.atom.count);
    return { mol.atom.position.x[i], mol.atom.position.y[i], mol.atom.position.z[i] };
}

inline i64 get_residue_atom_count(const MoleculeStructure& mol, ResIdx i) {
    ASSERT(0 <= i && i < mol.residue.count);
    return mol.residue.atom_range[i].ext();
}

inline soa_vec3 get_residue_positions(MoleculeStructure& mol, ResIdx i) {
    ASSERT(0 <= i && i < mol.residue.count);
    return {mol.atom.position.x + mol.residue.atom_range->beg, mol.atom.position.y + mol.residue.atom_range->beg,
            mol.atom.position.z + mol.residue.atom_range->beg};
}

inline const soa_vec3 get_residue_positions(const MoleculeStructure& mol, ResIdx i) {
    ASSERT(0 <= i && i < mol.residue.count);
    return {mol.atom.position.x + mol.residue.atom_range->beg, mol.atom.position.y + mol.residue.atom_range->beg,
            mol.atom.position.z + mol.residue.atom_range->beg};
}

inline const Element* get_residue_elements(const MoleculeStructure& mol, ResIdx i) {
    ASSERT(0 <= i && i < mol.residue.count);
    return mol.atom.element + mol.residue.atom_range->beg;
}

/*
// General accessors
inline Array<Residue> get_residues(MoleculeStructure& mol) { return mol.residues; }
inline Array<const Residue> get_residues(const MoleculeStructure& mol) { return mol.residues; }
inline Array<Chain> get_chains(MoleculeStructure& mol) { return mol.chains; }
inline Array<const Chain> get_chains(const MoleculeStructure& mol) { return mol.chains; }
inline Array<Bond> get_covalent_bonds(MoleculeStructure& mol) { return mol.covalent_bonds; }
inline Array<const Bond> get_covalent_bonds(const MoleculeStructure& mol) { return mol.covalent_bonds; }

// Single atom access
// @PERFORMANCE: Avoid this function
inline vec3 get_position_xyz(const MoleculeStructure& mol, AtomIdx idx) {
    ASSERT(0 <= idx && idx < mol.atom.count);
    return {mol.atom.position.x[idx], mol.atom.position.y[idx], mol.atom.position.z[idx]};
}

inline Array<float> get_positions_x(MoleculeStructure& mol) { return Array<float>(mol.atom.position.x, mol.atom.count); }
inline Array<const float> get_positions_x(const MoleculeStructure& mol) { return Array<const float>(mol.atom.position.x, mol.atom.count); }
inline Array<float> get_positions_y(MoleculeStructure& mol) { return Array<float>(mol.atom.position.y, mol.atom.count); }
inline Array<const float> get_positions_y(const MoleculeStructure& mol) { return Array<const float>(mol.atom.position.y, mol.atom.count); }
inline Array<float> get_positions_z(MoleculeStructure& mol) { return Array<float>(mol.atom.position.z, mol.atom.count); }
inline Array<const float> get_positions_z(const MoleculeStructure& mol) { return Array<const float>(mol.atom.position.z, mol.atom.count); }

inline Array<float> get_radii(MoleculeStructure& mol) { return Array<float>(mol.atom.radius, mol.atom.count); }

inline Array<Element> get_elements(MoleculeStructure& mol) { return Array<Element>(mol.atom.element, mol.atom.count); }
inline Array<const Element> get_elements(const MoleculeStructure& mol) { return Array<const Element>(mol.atom.element, mol.atom.count); }
inline Array<Label> get_labels(MoleculeStructure& mol) { return Array<Label>(mol.atom.label, mol.atom.count); }
inline Array<const Label> get_labels(const MoleculeStructure& mol) { return Array<const Label>(mol.atom.label, mol.atom.count); }
inline Array<ResIdx> get_residue_indices(MoleculeStructure& mol) { return Array<ResIdx>(mol.atom.res_idx, mol.atom.count); }
inline Array<const ResIdx> get_residue_indices(const MoleculeStructure& mol) { return Array<const ResIdx>(mol.atom.res_idx, mol.atom.count); }

// Backbone accessors
inline Array<BackboneSegment> get_backbone(MoleculeStructure& mol, BackboneSequence seq) {
    ASSERT(0 <= seq.beg && seq.end <= mol.backbone.segments.count);
    return mol.backbone.segments.subarray(seq);
}

inline Array<const BackboneSegment> get_backbone(const MoleculeStructure& mol, BackboneSequence seq) {
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
*/

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

/*
inline Array<Residue> get_residues(MoleculeStructure& mol, const Chain& chain) { return mol.residues.subarray(chain.res_range); }
inline Array<const Residue> get_residues(const MoleculeStructure& mol, const Chain& chain) { return mol.residues.subarray(chain.res_range); }

inline Array<Element> get_elements(MoleculeStructure& mol, Chain& chain) { return get_elements(mol).subarray(chain.atom_range); }
inline Array<const Element> get_elements(const MoleculeStructure& mol, const Chain& chain) { return get_elements(mol).subarray(chain.atom_range); }

inline Array<Label> get_labels(MoleculeStructure& mol, Chain& chain) { return get_labels(mol).subarray(chain.atom_range); }
inline Array<const Label> get_labels(const MoleculeStructure& mol, const Chain& chain) { return get_labels(mol).subarray(chain.atom_range); }

// Res func
inline Array<Element> get_elements(MoleculeStructure& mol, const Residue& res) { return get_elements(mol).subarray(res.atom_range); }
inline Array<const Element> get_elements(const MoleculeStructure& mol, const Residue& res) { return get_elements(mol).subarray(res.atom_range); }
inline Array<Label> get_labels(MoleculeStructure& mol, const Residue& res) { return get_labels(mol).subarray(res.atom_range); }
inline Array<const Label> get_labels(const MoleculeStructure& mol, const Residue& res) { return get_labels(mol).subarray(res.atom_range); }
inline Array<Bond> get_bonds(MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg);
}
inline Array<const Bond> get_bonds(const MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg);
}
inline Array<Bond> get_internal_bonds(MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg);
}
inline Array<const Bond> get_internal_bonds(const MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.subarray(res.bond_idx.beg_internal, res.bond_idx.end_internal - res.bond_idx.beg_internal);
}
*/

struct AtomDescriptor {
    float x, y, z;
    ResIdx residue_index;
    Label label;
    Element element;
};

struct ResidueDescriptor {
    Label name;
    AtomRange atom_range;
    ResIdx id;
};

struct ChainDescriptor {
    Label id;
    ResRange residue_range;
};

struct SecondaryStructureDescriptor {
    SecondaryStructure type;
    ResRange residue_range;
};

using BondDescriptor = Bond;

struct MoleculeStructureDescriptor {
    i32 num_atoms = 0;
    AtomDescriptor* atoms = nullptr;

    i32 num_residues = 0;
    ResidueDescriptor* residues = nullptr;

    i32 num_chains = 0;
    ChainDescriptor* chains = nullptr;

    i32 num_bonds = 0;
    BondDescriptor* bonds = nullptr;

    i32 num_secondary_structures = 0;
    SecondaryStructureDescriptor* secondary_structures = nullptr;
};

bool init_molecule_structure(MoleculeStructure* mol, const MoleculeStructureDescriptor& desc);
void free_molecule_structure(MoleculeStructure* mol);
