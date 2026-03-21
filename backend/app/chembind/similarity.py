# app/chembind/similarity.py
"""Tanimoto similarity search over user analysis history."""
from __future__ import annotations

from typing import Any, Dict, List

from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect

from .rdkit_safe import compute_morgan_fp, smiles_to_mol, RdkitLimits


def _bitstring_to_fp(bitstring: str) -> ExplicitBitVect:
    """Reconstruct an ExplicitBitVect from a '010110...' string."""
    n = len(bitstring)
    fp = ExplicitBitVect(n)
    for i, ch in enumerate(bitstring):
        if ch == "1":
            fp.SetBit(i)
    return fp


def tanimoto_search(
    query_smiles: str,
    history_docs: List[Dict[str, Any]],
    top_k: int = 10,
    min_sim: float = 0.3,
) -> List[Dict[str, Any]]:
    """
    Compute Tanimoto similarity between query and each doc that has a
    morganFp2048 field.  Return top_k results above min_sim, sorted desc.
    """
    limits = RdkitLimits()
    query_mol = smiles_to_mol(query_smiles, limits)
    query_fp = _bitstring_to_fp(compute_morgan_fp(query_mol))

    results: list[dict] = []
    for doc in history_docs:
        fp_str = doc.get("morganFp2048")
        if not fp_str:
            continue

        doc_fp = _bitstring_to_fp(fp_str)
        sim = DataStructs.TanimotoSimilarity(query_fp, doc_fp)

        if sim >= min_sim:
            results.append({
                "doc_id": doc.get("id", ""),
                "smiles": doc.get("smiles", ""),
                "similarity": round(float(sim), 4),
                "descriptors": doc.get("descriptors", {}),
            })

    results.sort(key=lambda r: r["similarity"], reverse=True)
    return results[:top_k]
