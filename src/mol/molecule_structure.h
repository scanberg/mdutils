#pragma once

#include <core/types.h>
#include <core/vector_types.h>
#include <core/string_utils.h>
#include <mol/element.h>

using AtomIdx   = i32;
using ResIdx    = i32;
using ChainIdx  = i32;
using SegIdx    = i32;
using BondIdx   = i32;
using AtomFlags = u8;

using AtomRange = Range<AtomIdx>;
using ResRange  = Range<ResIdx>;
using SegRange  = Range<SegIdx>;
using BondRange = Range<BondIdx>;

using HydrogenBondAcceptor = AtomIdx;

constexpr ChainIdx INVALID_CHAIN_IDX = -1;

enum class SecondaryStructure : u32 { Undefined = 0, Coil = 0x000000FF, Helix = 0x0000FF00, Sheet = 0x00FF0000 };

struct Bond {
    AtomIdx idx[2] = {0, 0};
};

struct BackboneAtoms {
    AtomIdx n_idx = 0;
    AtomIdx ca_idx = 0;
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

// Interface to access molecular data
struct MoleculeStructure {
    // SOA Layout for atom data
    struct {
        i64 count = 0;
        // Aligned data for vectorized operations
        soa_vec3 position = {nullptr, nullptr, nullptr};
        float* radius = nullptr;
        float* mass = nullptr;

        // Non aligned data
        Element* element = nullptr;
        const char** name = nullptr;
        ResIdx* res_idx = nullptr;
        ChainIdx* chain_idx = nullptr;
        uint8_t* flags = nullptr;
        //SeqIdx* seq_idx = nullptr;
    } atom;

    struct {
        i64 count = 0;
        Bond* bond = nullptr;
    } covalent_bond;

    // SOA layout for residue data
    struct {
        i64 count = 0;
        ResIdx* id = nullptr;
        const char** name = nullptr;
        AtomRange* atom_range = nullptr;

        struct {
            // External bonds are implicitly given by the range in complete but not in internal
            // Bonds to previous residue in chain is given by: [complete.beg, intra.beg]
            // Bonds to next residue in chain is given by: [intra.end, complete.end]
            BondRange* complete = nullptr;
            BondRange* intra = nullptr;
        } bond;

        struct { // This is an optional structure for amino-acid/nucleic backbones if the dataset contains such things
            BackboneAtoms* atoms = nullptr;
            BackboneAngle* angle = nullptr;
            SecondaryStructure* secondary_structure = nullptr;
        } backbone;
    } residue;

    // Chains represent connected residues (through covalent bonds) and are usually defined as a single character starting from 'A'
    struct {
        i64 count = 0;
        const char** id = nullptr;
        AtomRange* atom_range = nullptr;
        ResRange* residue_range = nullptr;
    } chain;

    /*
    struct {
        struct {
            i64 count = 0;
            //BackboneAtoms* segment = nullptr;
            BackboneAngle* angle = nullptr;
            //SecondaryStructure* secondary_structure = nullptr;
        } segment;

        struct {
            i64 count = 0;
            SegRange* segment_range = nullptr;
        } sequence;
    } backbone;
    */

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

    struct {
        // This is for storing interned strings
        void* strpool = nullptr;
    } internal;

    operator bool() const { return atom.count > 0; }
};

inline vec3 get_atom_position(const MoleculeStructure& mol, AtomIdx i) {
    (void)i;
    ASSERT(0 <= i && i <= mol.atom.count);
    return { mol.atom.position.x[i], mol.atom.position.y[i], mol.atom.position.z[i] };
}

inline Array<float>       get_atom_positions_x(MoleculeStructure& mol) { return {mol.atom.position.x, mol.atom.count}; }
inline Array<const float> get_atom_positions_x(const MoleculeStructure& mol) { return {mol.atom.position.x, mol.atom.count}; }

inline Array<float>       get_atom_positions_y(MoleculeStructure& mol) { return {mol.atom.position.y, mol.atom.count}; }
inline Array<const float> get_atom_positions_y(const MoleculeStructure& mol) { return {mol.atom.position.y, mol.atom.count}; }

inline Array<float>       get_atom_positions_z(MoleculeStructure& mol) { return {mol.atom.position.z, mol.atom.count}; }
inline Array<const float> get_atom_positions_z(const MoleculeStructure& mol) { return {mol.atom.position.z, mol.atom.count}; }

inline Array<float> get_radii(MoleculeStructure& mol) { return Array<float>(mol.atom.radius, mol.atom.count); }
inline Array<Element> get_elements(MoleculeStructure& mol) { return Array<Element>(mol.atom.element, mol.atom.count); }
inline Array<const Element> get_elements(const MoleculeStructure& mol) { return Array<const Element>(mol.atom.element, mol.atom.count); }
inline Array<const char*> get_names(MoleculeStructure& mol) { return Array<const char*>(mol.atom.name, mol.atom.count); }
inline Array<const char*> get_names(const MoleculeStructure& mol) { return Array<const char*>(mol.atom.name, mol.atom.count); }
inline Array<ResIdx> get_residue_indices(MoleculeStructure& mol) { return Array<ResIdx>(mol.atom.res_idx, mol.atom.count); }
inline Array<const ResIdx> get_residue_indices(const MoleculeStructure& mol) { return Array<const ResIdx>(mol.atom.res_idx, mol.atom.count); }

inline Array<const char*> get_residue_names(const MoleculeStructure& mol) {
    return {mol.residue.name, mol.residue.count};
}

inline Array<ResIdx> get_reisidue_ids(MoleculeStructure& mol) {
    return {mol.residue.id, mol.residue.count};
}

inline Array<const ResIdx> get_reisidue_ids(const MoleculeStructure& mol) {
    return {mol.residue.id, mol.residue.count};
}

inline i64 get_residue_atom_count(const MoleculeStructure& mol, ResIdx i) {
    (void)i;
    ASSERT(0 <= i && i < mol.residue.count);
    return mol.residue.atom_range[i].ext();
}

inline soa_vec3 get_residue_positions(MoleculeStructure& mol, ResIdx i) {
    (void)i;
    ASSERT(0 <= i && i < mol.residue.count);
    return mol.atom.position + mol.residue.atom_range->beg;
}

inline const soa_vec3 get_residue_positions(const MoleculeStructure& mol, ResIdx i) {
    (void)i;
    ASSERT(0 <= i && i < mol.residue.count);
    return mol.atom.position + mol.residue.atom_range->beg;
}

inline Array<Element> get_residue_elements(MoleculeStructure& mol, ResIdx i) {
    ASSERT(0 <= i && i < mol.residue.count);
    return {mol.atom.element + mol.residue.atom_range[i].beg, mol.residue.atom_range[i].ext()};
}

inline Array<const Element> get_residue_elements(const MoleculeStructure& mol, ResIdx i) {
    ASSERT(0 <= i && i < mol.residue.count);
    return {mol.atom.element + mol.residue.atom_range[i].beg, mol.residue.atom_range[i].ext()};
}

inline Array<const ResRange> get_chain_residue_ranges(const MoleculeStructure& mol) {
    return {mol.chain.residue_range, mol.chain.count};
}

inline Array<const BackboneAtoms> get_residue_backbone_atoms(const MoleculeStructure& mol) {
    if (mol.residue.backbone.atoms) {
        return {mol.residue.backbone.atoms, mol.residue.count};
    }
    return {};
}

inline Array<const BackboneAtoms> get_residue_backbone_atoms(const MoleculeStructure& mol, ChainIdx i) {
    ASSERT(0 <= i && i < mol.chain.count);
    if (mol.residue.backbone.atoms) {
        const auto range = mol.chain.residue_range[i];
        return {mol.residue.backbone.atoms + range.beg, range.ext()};
    }
    return {};
}

inline Array<const BackboneAngle> get_residue_backbone_angles(const MoleculeStructure& mol) {
    if (mol.residue.backbone.angle) {
        return {mol.residue.backbone.angle, mol.residue.count};
    }
    return {};
}

inline Array<const BackboneAngle> get_residue_backbone_angles(const MoleculeStructure& mol, ChainIdx i) {
    ASSERT(0 <= i && i < mol.chain.count);
    if (mol.residue.backbone.atoms) {
        const auto range = mol.chain.residue_range[i];
        return {mol.residue.backbone.angle + range.beg, range.ext()};
    }
    return {};
}


inline Array<BackboneAngle> get_residue_backbone_angles(MoleculeStructure& mol, ChainIdx i) {
    ASSERT(0 <= i && i < mol.chain.count);
    if (mol.residue.backbone.atoms) {
        const auto range = mol.chain.residue_range[i];
        return {mol.residue.backbone.angle + range.beg, range.ext()};
    }
    return {};
}

/*
inline Array<BackboneAtoms> get_backbone_segments(MoleculeStructure& mol, ChainIdx i) {
    return {mol.residue.backbone.atoms + mol.chain.residue_range[i].beg, mol.chain.residue_range[i].ext()};
}

inline Array<const BackboneAtoms> get_backbone_segments(const MoleculeStructure& mol, ChainIdx i) {
    return {mol.backbone.segment.segment + seq.beg, seq.ext()};
}

inline Array<BackboneAngle> get_backbone_angles(MoleculeStructure& mol, BackboneSequence seq) {
    return {mol.backbone.segment.angle + seq.beg, seq.ext()};
}

inline const Array<BackboneSequence> get_backbone_sequences(const MoleculeStructure& mol) {
    return {mol.backbone.sequence.segment_range, mol.backbone.sequence.count};
}
*/

struct AtomDescriptor {
    float x, y, z;
    ResIdx residue_index;
    CStringView name;
    Element element;
};

struct ResidueDescriptor {
    CStringView name;
    AtomRange atom_range;
    ResIdx id;
};

struct ChainDescriptor {
    CStringView id;
    ResRange residue_range;
};

struct SecondaryStructureDescriptor {
    SecondaryStructure type;
    ResRange residue_range;
};

using BondDescriptor = Bond;

struct MoleculeStructureDescriptor {
    i64 num_atoms = 0;
    AtomDescriptor* atoms = nullptr;

    i64 num_residues = 0;
    ResidueDescriptor* residues = nullptr;

    i64 num_chains = 0;
    ChainDescriptor* chains = nullptr;

    i64 num_bonds = 0;
    BondDescriptor* bonds = nullptr;

    i64 num_secondary_structures = 0;
    SecondaryStructureDescriptor* secondary_structures = nullptr;
};

bool init_molecule_structure(MoleculeStructure* mol, const MoleculeStructureDescriptor& desc);
void free_molecule_structure(MoleculeStructure* mol);
