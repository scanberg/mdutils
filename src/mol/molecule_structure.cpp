#include "molecule_structure.h"

// 32-byte alignment for 256-bit vectorization (AVX+ architectures)
#define ALIGNMENT 32

bool init_molecule_structure(MoleculeStructure* mol, int32 num_atoms, int32 num_bonds, int32 num_residues, int32 num_chains, int32 num_sequences, int32 num_backbone_segments, int32 num_backbone_sequences,
                             int32 num_hydrogen_bond_donors, int32 num_hydrogen_bond_acceptors) {
    free_molecule_structure(mol);

    int64 alloc_size = 0;
    alloc_size += num_atoms * (sizeof(Element) + sizeof(Label) + sizeof(ResIdx) + sizeof(ChainIdx) + sizeof(SeqIdx));
    alloc_size += num_bonds * sizeof(Bond);
    alloc_size += num_residues * sizeof(Residue);
    alloc_size += num_chains * sizeof(Chain);
    alloc_size += num_sequences * sizeof(Sequence);
    alloc_size += num_backbone_segments * sizeof(BackboneSegment);
    alloc_size += num_backbone_segments * sizeof(BackboneAngle);
    alloc_size += num_backbone_sequences * sizeof(BackboneSequence);
    alloc_size += num_hydrogen_bond_donors * sizeof(HydrogenBondDonor);
    alloc_size += num_hydrogen_bond_acceptors * sizeof(HydrogenBondAcceptor);

    // Allocate Aligned data (@NOTE: Is perhaps not necessary as trajectory data is not aligned anyways...)
    void* data = MALLOC(alloc_size);
    void* pos_data = ALIGNED_MALLOC((num_atoms * sizeof(float) + ALIGNMENT) * 3, ALIGNMENT);
    void* vel_data = ALIGNED_MALLOC((num_atoms * sizeof(float) + ALIGNMENT) * 3, ALIGNMENT);
    void* rad_data = ALIGNED_MALLOC((num_atoms * sizeof(float) + ALIGNMENT) * 1, ALIGNMENT);
    void* mas_data = ALIGNED_MALLOC((num_atoms * sizeof(float) + ALIGNMENT) * 1, ALIGNMENT);

    if (!data) return false;
    if (!pos_data) return false;
    if (!vel_data) return false;
    if (!rad_data) return false;
    if (!mas_data) return false;

    memset(data, 0, alloc_size);
    memset(pos_data, 0, (num_atoms * sizeof(float) + ALIGNMENT) * 3);
    memset(vel_data, 0, (num_atoms * sizeof(float) + ALIGNMENT) * 3);
    memset(rad_data, 0, (num_atoms * sizeof(float) + ALIGNMENT) * 1);
    memset(mas_data, 0, (num_atoms * sizeof(float) + ALIGNMENT) * 1);

    mol->atom.count = num_atoms;

    mol->atom.position.x = (float*)pos_data;
    mol->atom.position.y = (float*)get_next_aligned_adress(mol->atom.position.x + num_atoms, ALIGNMENT);
    mol->atom.position.z = (float*)get_next_aligned_adress(mol->atom.position.y + num_atoms, ALIGNMENT);
    ASSERT(IS_ALIGNED(mol->atom.position.x, ALIGNMENT));
    ASSERT(IS_ALIGNED(mol->atom.position.y, ALIGNMENT));
    ASSERT(IS_ALIGNED(mol->atom.position.z, ALIGNMENT));

    mol->atom.velocity.x = (float*)vel_data;
    mol->atom.velocity.y = (float*)get_next_aligned_adress(mol->atom.velocity.x + num_atoms, ALIGNMENT);
    mol->atom.velocity.z = (float*)get_next_aligned_adress(mol->atom.velocity.y + num_atoms, ALIGNMENT);
    ASSERT(IS_ALIGNED(mol->atom.velocity.x, ALIGNMENT));
    ASSERT(IS_ALIGNED(mol->atom.velocity.y, ALIGNMENT));
    ASSERT(IS_ALIGNED(mol->atom.velocity.z, ALIGNMENT));

    mol->atom.radius = (float*)rad_data;
    ASSERT(IS_ALIGNED(mol->atom.radius, ALIGNMENT));

    mol->atom.mass = (float*)mas_data;
    ASSERT(IS_ALIGNED(mol->atom.mass, ALIGNMENT));

    mol->atom.element = (Element*)data;
    mol->atom.label = (Label*)(mol->atom.element + num_atoms);
    mol->atom.res_idx = (ResIdx*)(mol->atom.label + num_atoms);
    mol->atom.chain_idx = (ChainIdx*)(mol->atom.res_idx + num_atoms);
    mol->atom.seq_idx = (ChainIdx*)(mol->atom.chain_idx + num_atoms);

    mol->covalent_bonds = {(Bond*)(mol->atom.seq_idx + num_atoms), num_bonds};
    mol->residues = {(Residue*)(mol->covalent_bonds.end()), num_residues};
    mol->chains = {(Chain*)(mol->residues.end()), num_chains};
    mol->sequences = {(Sequence*)(mol->chains.end()), num_sequences};
    mol->backbone.segments = {(BackboneSegment*)(mol->sequences.end()), num_backbone_segments};
    mol->backbone.angles = {(BackboneAngle*)(mol->backbone.segments.end()), num_backbone_segments};
    mol->backbone.sequences = {(BackboneSequence*)(mol->backbone.angles.end()), num_backbone_sequences};
    mol->hydrogen_bond.donors = {(HydrogenBondDonor*)(mol->backbone.sequences.end()), num_hydrogen_bond_donors};
    mol->hydrogen_bond.acceptors = {(HydrogenBondAcceptor*)(mol->hydrogen_bond.donors.end()), num_hydrogen_bond_acceptors};

    return true;
}

void free_molecule_structure(MoleculeStructure* mol) {
    ASSERT(mol);
    if (mol->atom.position.x) ALIGNED_FREE(mol->atom.position.x);
    if (mol->atom.velocity.x) ALIGNED_FREE(mol->atom.velocity.x);
    if (mol->atom.radius) ALIGNED_FREE(mol->atom.radius);
    if (mol->atom.mass) ALIGNED_FREE(mol->atom.mass);
    if (mol->atom.element) FREE(mol->atom.element);
    *mol = {};
}
