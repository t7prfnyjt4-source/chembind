# app/chembind/conformers.py
"""Conformer ensemble generation using RDKit ETKDGv3 + MMFF optimization."""
from __future__ import annotations

import hashlib
import json
from typing import Any, Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from .rdkit_safe import smiles_to_mol, SmilesValidationError, RdkitLimits


def generate_conformers(
    smiles: str,
    num_confs: int = 20,
    max_atoms: int = 80,
    seed: int = 42,
) -> List[Dict[str, Any]]:
    """
    Generate a conformer ensemble for a molecule.

    Returns list of dicts sorted by energy (ascending):
      [{conf_id, energy, atoms: [{symbol, x, y, z}]}]
    """
    limits = RdkitLimits(max_atoms=max_atoms)
    mol = smiles_to_mol(smiles, limits)

    mol_h = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1

    cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=params)
    if not cids:
        raise SmilesValidationError("Could not generate 3D conformers for this molecule")

    # MMFF optimize each conformer
    results_raw = AllChem.MMFFOptimizeMoleculeConfs(mol_h, numThreads=1)

    conformers: List[Dict[str, Any]] = []
    for idx, cid in enumerate(cids):
        converged, energy = results_raw[idx] if idx < len(results_raw) else (1, 0.0)

        conf = mol_h.GetConformer(cid)
        atoms = []
        for i in range(mol_h.GetNumAtoms()):
            atom = mol_h.GetAtomWithIdx(i)
            # Skip hydrogens in output to keep payload small
            if atom.GetAtomicNum() == 1:
                continue
            pos = conf.GetAtomPosition(i)
            atoms.append({
                "symbol": atom.GetSymbol(),
                "x": round(pos.x, 4),
                "y": round(pos.y, 4),
                "z": round(pos.z, 4),
            })

        conformers.append({
            "conf_id": int(cid),
            "energy": round(float(energy), 4),
            "atoms": atoms,
        })

    conformers.sort(key=lambda c: c["energy"])
    return conformers


def conformer_cache_key(smiles: str, num_confs: int) -> str:
    """Deterministic cache key for conformer results."""
    raw = f"{smiles}:{num_confs}"
    return f"conf:{hashlib.sha256(raw.encode()).hexdigest()[:32]}"
