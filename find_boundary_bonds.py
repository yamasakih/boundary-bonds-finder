from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol


def find_boundary_bonds(mol, scaffold_atom_indices=None):
    if not scaffold_atom_indices:
        scaffold = GetScaffoldForMol(mol)    
        scaffold_atom_indices = mol.GetSubstructMatch(scaffold)        
    return [bond for atom_idx in scaffold_atom_indices for bond in mol.GetAtomWithIdx(atom_idx).GetBonds()
                                                                   if bond.GetOtherAtomIdx(atom_idx) not in scaffold_atom_indices]
