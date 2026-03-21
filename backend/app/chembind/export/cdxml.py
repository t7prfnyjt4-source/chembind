# app/chembind/export/cdxml.py
"""Generate CDXML (ChemDraw XML) from RDKit Mol objects."""
from __future__ import annotations

import xml.etree.ElementTree as ET

from rdkit import Chem
from rdkit.Chem import AllChem


# Scale factor: RDKit coords (Angstroms) → CDXML points
SCALE = 30.0
OFFSET_X = 200.0
OFFSET_Y = 200.0

BOND_ORDER_MAP = {
    Chem.BondType.SINGLE: "1",
    Chem.BondType.DOUBLE: "2",
    Chem.BondType.TRIPLE: "3",
    Chem.BondType.AROMATIC: "1.5",
}


def mol_to_cdxml(mol: Chem.Mol) -> str:
    """
    Convert an RDKit Mol to CDXML string.
    Structure-only (no reaction schemes, R-groups).
    """
    # Compute 2D coordinates
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()

    # Build XML tree
    cdxml = ET.Element("CDXML")
    page = ET.SubElement(cdxml, "page")
    fragment = ET.SubElement(page, "fragment")

    # Atoms
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()

        n = ET.SubElement(fragment, "n", {
            "id": str(i + 1),
            "p": f"{pos.x * SCALE + OFFSET_X:.4f} {-pos.y * SCALE + OFFSET_Y:.4f}",
            "Element": str(atom.GetAtomicNum()),
        })

        # Non-carbon atoms get text labels
        if symbol != "C":
            t = ET.SubElement(n, "t")
            s = ET.SubElement(t, "s")
            s.text = symbol

    # Bonds
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx() + 1
        end = bond.GetEndAtomIdx() + 1
        order = BOND_ORDER_MAP.get(bond.GetBondType(), "1")

        ET.SubElement(fragment, "b", {
            "B": str(begin),
            "E": str(end),
            "Order": order,
        })

    # Serialize
    return ET.tostring(cdxml, encoding="unicode", xml_declaration=True)
