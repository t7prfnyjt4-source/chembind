# app/chembind/docking/preparation.py
"""Ligand and receptor preparation for docking."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_ligand_from_smiles(smiles: str, tmpdir: str) -> str:
    """
    Convert SMILES → 3D SDF → PDBQT via Meeko.
    Returns path to ligand PDBQT file.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol_h = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol_h, params)
    if status != 0:
        raise ValueError("Could not generate 3D coordinates for ligand")

    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=200)

    sdf_path = os.path.join(tmpdir, "ligand.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol_h)
    writer.close()

    # Meeko: SDF → PDBQT
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol_h)

    pdbqt_path = os.path.join(tmpdir, "ligand.pdbqt")
    for setup in mol_setups:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        if is_ok:
            with open(pdbqt_path, "w") as f:
                f.write(pdbqt_string)
            return pdbqt_path

    raise ValueError("Meeko failed to produce ligand PDBQT")


def prepare_receptor_from_pdb(pdb_content: str, tmpdir: str) -> str:
    """
    Write PDB content to file and convert to PDBQT.
    Returns path to receptor PDBQT file.

    Uses a simple approach: strip HETATM lines and convert via
    basic PDB → PDBQT formatting. For production, use ADFR suite's
    mk_prepare_receptor.py if available.
    """
    pdb_path = os.path.join(tmpdir, "receptor.pdb")
    with open(pdb_path, "w") as f:
        f.write(pdb_content)

    pdbqt_path = os.path.join(tmpdir, "receptor.pdbqt")

    # Try mk_prepare_receptor.py first (ADFR suite)
    try:
        result = subprocess.run(
            ["mk_prepare_receptor", "-i", pdb_path, "-o", pdbqt_path],
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode == 0 and os.path.exists(pdbqt_path):
            return pdbqt_path
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    # Fallback: simple PDB → PDBQT conversion
    # Keep only ATOM lines, add dummy charges
    lines = []
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("TER") or line.startswith("END"):
            if line.startswith("ATOM"):
                # Pad to 80 chars, add autodock type
                padded = line.ljust(77)
                element = line[76:78].strip() if len(line) > 76 else line[12:14].strip()
                element = element[0] if element else "C"
                pdbqt_line = f"{padded[:54]}  0.00  0.00    +0.000 {element:>2s}"
                lines.append(pdbqt_line)
            else:
                lines.append(line)

    if not lines:
        raise ValueError("PDB file contains no ATOM records")

    with open(pdbqt_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return pdbqt_path
