# app/chembind/export/formats.py
"""Export molecules to MOL and SDF formats."""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem


def mol_to_mol(mol: Chem.Mol) -> str:
    """Export as 2D MOL block (V2000)."""
    AllChem.Compute2DCoords(mol)
    return Chem.MolToMolBlock(mol)


def mol_to_sdf(mol: Chem.Mol) -> str:
    """Export as 3D SDF block (with MMFF-optimized coords)."""
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol_h, params)
    if status != 0:
        # Fallback to 2D if 3D embedding fails
        AllChem.Compute2DCoords(mol_h)
    else:
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=200)

    # Remove Hs for cleaner output
    mol_out = Chem.RemoveHs(mol_h)
    block = Chem.MolToMolBlock(mol_out)
    return block + "\n$$$$\n"
