# app/chembind/docking/engine.py
"""AutoDock Vina docking engine wrapper."""
from __future__ import annotations

import os
import re
from typing import Any, Dict, List


def parse_pdbqt_poses(pdbqt_path: str) -> List[Dict[str, Any]]:
    """
    Parse multi-model PDBQT output into list of poses.
    Each pose: {pose_id, score, atoms: [{name, x, y, z, element}]}
    """
    poses: List[Dict[str, Any]] = []
    current_atoms: List[Dict[str, Any]] = []
    current_score = 0.0
    pose_id = 0

    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                pose_id += 1
                current_atoms = []
                current_score = 0.0
            elif line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        current_score = float(parts[3])
                    except ValueError:
                        pass
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip() if len(line) > 76 else name[0]
                current_atoms.append({
                    "name": name,
                    "x": round(x, 4),
                    "y": round(y, 4),
                    "z": round(z, 4),
                    "element": element,
                })
            elif line.startswith("ENDMDL"):
                poses.append({
                    "pose_id": pose_id,
                    "score": current_score,
                    "atoms": current_atoms,
                })

    # If no MODEL/ENDMDL, treat as single pose
    if not poses and current_atoms:
        poses.append({
            "pose_id": 1,
            "score": current_score,
            "atoms": current_atoms,
        })

    return poses


def run_vina_docking(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    center: tuple[float, float, float],
    size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    num_modes: int = 9,
    exhaustiveness: int = 8,
    output_dir: str = ".",
) -> List[Dict[str, Any]]:
    """
    Run AutoDock Vina docking and return parsed poses.
    """
    from vina import Vina

    v = Vina(sf_name="vina")
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)

    v.compute_vina_maps(
        center=list(center),
        box_size=list(size),
    )

    v.dock(
        exhaustiveness=exhaustiveness,
        n_poses=num_modes,
    )

    output_path = os.path.join(output_dir, "poses.pdbqt")
    v.write_poses(output_path, n_poses=num_modes, overwrite=True)

    return parse_pdbqt_poses(output_path)
