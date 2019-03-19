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

using BackboneSequence = ResRange;
using BackboneAngle = vec2;

struct HydrogenBondDonor {
    AtomIdx donor_idx = 0;
    AtomIdx hydro_idx = 0;
};

using HydrogenBondAcceptor = AtomIdx;

struct Residue {
    Label name{};
    ResIdx id = -1;
    ChainIdx chain_idx = -1;

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

struct Chain {
    Label id{};

    ResRange res_range;
    AtomRange atom_range;
};

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
    } atom;

    Array<Bond> covalent_bonds{};
    Array<Residue> residues{};
    Array<Chain> chains{};

    struct {
        // Segments and angles should match in length and if not zero, they should match the number of residues
        Array<BackboneSegment> segments{};
        Array<BackboneAngle> angles{};
        Array<BackboneSequence> sequences{};
    } backbone;

    struct {
        Array<HydrogenBondDonor> donors{};
        Array<HydrogenBondAcceptor> acceptors{};
    } hydrogen_bond;

    operator bool() const { return atom.count > 0; }
};

// General accessors
inline Array<Residue> get_residues(MoleculeStructure& mol) { return mol.residues; }
inline Array<const Residue> get_residues(const MoleculeStructure& mol) { return mol.residues; }
inline Array<Chain> get_chains(MoleculeStructure& mol) { return mol.chains; }
inline Array<const Chain> get_chains(const MoleculeStructure& mol) { return mol.chains; }
inline Array<Bond> get_covalent_bonds(MoleculeStructure& mol) { return mol.covalent_bonds; }
inline Array<const Bond> get_covalent_bonds(const MoleculeStructure& mol) { return mol.covalent_bonds; }

inline Float3Stream get_position_stream(MoleculeStructure& mol) { return { mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, mol.atom.count }; }
inline Float3Stream get_velocity_stream(MoleculeStructure& mol) { return { mol.atom.velocity.x, mol.atom.velocity.y, mol.atom.velocity.z, mol.atom.count }; }

// Single atom access
inline vec3 get_position_xyz(const MoleculeStructure& mol, AtomIdx idx) {
	ASSERT(0 <= idx && idx < mol.atom.count);
	return { mol.atom.position.x[idx], mol.atom.position.y[idx], mol.atom.position.z[idx] };
}

inline Array<float> get_positions_x(MoleculeStructure& mol) { return Array<float>(mol.atom.position.x, mol.atom.count); }
inline Array<const float> get_positions_x(const MoleculeStructure& mol) { return Array<const float>(mol.atom.position.x, mol.atom.count); }
inline Array<float> get_positions_y(MoleculeStructure& mol) { return Array<float>(mol.atom.position.y, mol.atom.count); }
inline Array<const float> get_positions_y(const MoleculeStructure& mol) { return Array<const float>(mol.atom.position.y, mol.atom.count); }
inline Array<float> get_positions_z(MoleculeStructure& mol) { return Array<float>(mol.atom.position.z, mol.atom.count); }
inline Array<const float> get_positions_z(const MoleculeStructure& mol) { return Array<const float>(mol.atom.position.z, mol.atom.count); }

inline Array<float> get_velocities_x(MoleculeStructure& mol) { return Array<float>(mol.atom.velocity.x, mol.atom.count); }
inline Array<const float> get_velocities_x(const MoleculeStructure& mol) { return Array<const float>(mol.atom.velocity.x, mol.atom.count); }
inline Array<float> get_velocities_y(MoleculeStructure& mol) { return Array<float>(mol.atom.velocity.y, mol.atom.count); }
inline Array<const float> get_velocities_y(const MoleculeStructure& mol) { return Array<const float>(mol.atom.velocity.y, mol.atom.count); }
inline Array<float> get_velocities_z(MoleculeStructure& mol) { return Array<float>(mol.atom.velocity.z, mol.atom.count); }
inline Array<const float> get_velocities_z(const MoleculeStructure& mol) { return Array<const float>(mol.atom.velocity.z, mol.atom.count); }

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
    return mol.backbone.segments.subarray(seq.beg, seq.end - seq.beg);
}

inline Array<const BackboneSegment> get_backbone(const MoleculeStructure& mol, BackboneSequence seq) {
    ASSERT(0 <= seq.beg && seq.end <= mol.backbone.segments.count);
    return mol.backbone.segments.subarray(seq.beg, seq.end - seq.beg);
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

inline Array<Residue> get_residues(MoleculeStructure& mol, const Chain& chain) { return mol.residues.subarray(chain.res_range); }
inline Array<const Residue> get_residues(const MoleculeStructure& mol, const Chain& chain) { return mol.residues.subarray(chain.res_range); }

inline Array<Element> get_elements(MoleculeStructure& mol, Chain& chain) {
    return get_elements(mol).subarray(chain.atom_range);
}

inline Array<const Element> get_elements(const MoleculeStructure& mol, const Chain& chain) {
    return get_elements(mol).subarray(chain.atom_range);
}

inline Array<Label> get_labels(MoleculeStructure& mol, Chain& chain) {
    return get_labels(mol).subarray(chain.atom_range);
}

inline Array<const Label> get_labels(const MoleculeStructure& mol, const Chain& chain) {
    return get_labels(mol).subarray(chain.atom_range);
}

// Res func
inline Array<Element> get_elements(MoleculeStructure& mol, const Residue& res) { return get_elements(mol).subarray(res.atom_range); }
inline Array<const Element> get_elements(const MoleculeStructure& mol, const Residue& res) { return get_elements(mol).subarray(res.atom_range); }
inline Array<Label> get_labels(MoleculeStructure& mol, const Residue& res) { return get_labels(mol).subarray(res.atom_range); }
inline Array<const Label> get_labels(const MoleculeStructure& mol, const Residue& res) { return get_labels(mol).subarray(res.atom_range); }
inline Array<Bond> get_bonds(MoleculeStructure& mol, const Residue& res) { return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg); }
inline Array<const Bond> get_bonds(const MoleculeStructure& mol, const Residue& res) { return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg); }
inline Array<Bond> get_internal_bonds(MoleculeStructure& mol, const Residue& res) { return mol.covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg); }
inline Array<const Bond> get_internal_bonds(const MoleculeStructure& mol, const Residue& res) {
    return mol.covalent_bonds.subarray(res.bond_idx.beg_internal, res.bond_idx.end_internal - res.bond_idx.beg_internal);
}

bool init_molecule_structure(MoleculeStructure* mol, int32 num_atoms, int32 num_bonds, int32 num_residues, int32 num_chains, int32 num_backbone_segments = 0, int32 num_backbone_sequences = 0,
                             int32 num_hydrogen_bond_donors = 0, int32 num_hydrogen_bond_acceptors = 0);
void free_molecule_structure(MoleculeStructure* mol);
