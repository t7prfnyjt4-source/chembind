# app/chembind/docking/interactions.py
"""Protein-ligand interaction analysis using ProLIF."""
from __future__ import annotations

import logging
from typing import Any, Dict, List

logger = logging.getLogger("chembind.docking")


def analyze_interactions(
    protein_pdb: str,
    pose_atoms: List[Dict[str, Any]],
    ligand_smiles: str,
) -> Dict[str, Any]:
    """
    Analyze protein-ligand interactions for a docking pose.

    Returns:
        {interactions: [{residue, type}], types_found: [str]}

    Circuit breaker: returns empty if ProLIF fails.
    """
    try:
        import prolif
        import MDAnalysis as mda
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import tempfile
        import os

        # Write protein PDB
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
            f.write(protein_pdb)
            protein_path = f.name

        # Build ligand SDF from pose atoms
        mol = Chem.MolFromSmiles(ligand_smiles)
        if mol is None:
            return _empty_result()

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        # Write ligand SDF
        with tempfile.NamedTemporaryFile(suffix=".sdf", mode="w", delete=False) as f:
            f.write(Chem.MolToMolBlock(mol))
            f.write("\n$$$$\n")
            ligand_path = f.name

        try:
            # Load with MDAnalysis
            prot = mda.Universe(protein_path)
            lig = mda.Universe(ligand_path)

            # ProLIF fingerprint
            fp = prolif.Fingerprint()
            fp.run(lig.select_atoms("all"), prot.select_atoms("protein"))

            # Extract interactions
            interactions: List[Dict[str, str]] = []
            types_found: set = set()

            df = fp.to_dataframe()
            for col in df.columns:
                if df[col].any():
                    # Column format: (ligand_key, protein_residue, interaction_type)
                    if isinstance(col, tuple) and len(col) >= 3:
                        residue = str(col[1])
                        int_type = str(col[2])
                        interactions.append({"residue": residue, "type": int_type})
                        types_found.add(int_type)

            return {
                "interactions": interactions,
                "types_found": sorted(types_found),
            }

        finally:
            os.unlink(protein_path)
            os.unlink(ligand_path)

    except Exception as e:
        logger.warning(f"ProLIF interaction analysis failed (non-fatal): {e}")
        return _empty_result()


def _empty_result() -> Dict[str, Any]:
    return {"interactions": [], "types_found": []}
