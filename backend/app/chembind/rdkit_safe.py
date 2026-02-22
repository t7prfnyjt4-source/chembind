# app/chembind/rdkit_safe.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.error")

class SmilesValidationError(ValueError):
    pass


@dataclass(frozen=True)
class RdkitLimits:
    max_smiles_len: int = 500
    max_atoms: int = 200  # tune later; 200 is a reasonable MVP guard


def _ensure_limits(smiles: str, limits: RdkitLimits) -> None:
    if not isinstance(smiles, str) or not smiles.strip():
        raise SmilesValidationError("SMILES is empty")
    if len(smiles) > limits.max_smiles_len:
        raise SmilesValidationError(f"SMILES too long (>{limits.max_smiles_len})")


def smiles_to_mol(smiles: str, limits: RdkitLimits) -> Chem.Mol:
    _ensure_limits(smiles, limits)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesValidationError("Invalid SMILES")

    atom_count = mol.GetNumAtoms()
    if atom_count > limits.max_atoms:
        raise SmilesValidationError(f"Too many atoms ({atom_count} > {limits.max_atoms})")

    return mol


def compute_descriptors(smiles: str, limits: RdkitLimits | None = None) -> Dict[str, Any]:
    limits = limits or RdkitLimits()
    mol = smiles_to_mol(smiles, limits)

    # keep stable + predictable outputs
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)

    return {
        "mw": float(mw),
        "logp": float(logp),
        "hbd": int(hbd),
        "hba": int(hba),
        "tpsa": float(tpsa),
        "atom_count": int(mol.GetNumAtoms()),
        "heavy_atom_count": int(heavy_atoms),
        "ring_count": int(rings),
    }
