#include "molecule_structure.h"

#define ALIGNMENT 64

bool init_molecule_structure(MoleculeStructure* mol, int32 num_atoms, int32 num_bonds, int32 num_residues, int32 num_chains, int32 num_backbone_segments, int32 num_backbone_sequences,
                             int32 num_hydrogen_bond_donors, int32 num_hydrogen_bond_acceptors) {
    free_molecule_structure(mol);

    int64 alloc_size = 0;
    alloc_size += num_atoms * (sizeof(Element) + sizeof(Label) + sizeof(ResIdx));
    alloc_size += (num_atoms * sizeof(float) + ALIGNMENT) * (3 + 3 + 1);  // extra padding for alignment (pos, vel, radius)

    alloc_size += num_bonds * sizeof(Bond);
    alloc_size += num_residues * sizeof(Residue);
    alloc_size += num_chains * sizeof(Chain);
    alloc_size += num_backbone_segments * sizeof(BackboneSegment);
    alloc_size += num_backbone_segments * sizeof(BackboneAngle);
    alloc_size += num_backbone_sequences * sizeof(BackboneSequence);
    alloc_size += num_hydrogen_bond_donors * sizeof(HydrogenBondDonor);
    alloc_size += num_hydrogen_bond_acceptors * sizeof(HydrogenBondAcceptor);

    // Aligned

    void* data = MALLOC(alloc_size);
    if (!data) return false;

    mol->atom.count = num_atoms;

    mol->atom.element = (Element*)data;
    mol->atom.label = (Label*)get_elements(*mol).end();
    mol->atom.res_idx = (ResIdx*)get_labels(*mol).end();

    mol->atom.position.x = (float*)get_next_aligned_adress(mol->atom.res_idx + num_atoms, ALIGNMENT);
    mol->atom.position.y = (float*)get_next_aligned_adress(mol->atom.position.x + num_atoms, ALIGNMENT);
    mol->atom.position.z = (float*)get_next_aligned_adress(mol->atom.position.y + num_atoms, ALIGNMENT);

    mol->atom.velocity.x = (float*)get_next_aligned_adress(mol->atom.position.z + num_atoms, ALIGNMENT);
    mol->atom.velocity.y = (float*)get_next_aligned_adress(mol->atom.velocity.x + num_atoms, ALIGNMENT);
    mol->atom.velocity.z = (float*)get_next_aligned_adress(mol->atom.velocity.y + num_atoms, ALIGNMENT);

    mol->atom.radius = (float*)get_next_aligned_adress(mol->atom.velocity.z + num_atoms, ALIGNMENT);

    mol->covalent_bonds = {(Bond*)(mol->atom.radius + num_atoms), num_bonds};
    mol->residues = {(Residue*)(mol->covalent_bonds.end()), num_residues};
    mol->chains = {(Chain*)(mol->residues.end()), num_chains};
    mol->backbone.segments = {(BackboneSegment*)(mol->chains.end()), num_backbone_segments};
    mol->backbone.angles = {(BackboneAngle*)(mol->backbone.segments.end()), num_backbone_segments};
    mol->backbone.sequences = {(BackboneSequence*)(mol->backbone.angles.end()), num_backbone_sequences};
    mol->hydrogen_bond.donors = {(HydrogenBondDonor*)(mol->backbone.sequences.end()), num_hydrogen_bond_donors};
    mol->hydrogen_bond.acceptors = {(HydrogenBondAcceptor*)(mol->hydrogen_bond.donors.end()), num_hydrogen_bond_acceptors};

    return true;
}

void free_molecule_structure(MoleculeStructure* mol) {
    ASSERT(mol);
    if (mol->atom.element) FREE(mol->atom.element);
    *mol = {};
}
