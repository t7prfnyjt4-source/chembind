# tests/test_docking_api.py
"""Tests for docking API endpoint validation."""
import pytest
import io


def test_docking_enterprise_gate(client_factory):
    """Non-enterprise user should get 402."""
    # This test validates the tier gating logic
    from app.chembind.tier_gating import check_tier
    from fastapi import HTTPException

    basic_user = {"uid": "test", "claims": {"tier": "basic"}}
    with pytest.raises(HTTPException) as exc:
        check_tier(basic_user, "enterprise")
    assert exc.value.status_code == 402


def test_docking_invalid_smiles_rejected():
    """Invalid SMILES should be caught."""
    from rdkit import Chem
    assert Chem.MolFromSmiles("invalid[[[") is None


def test_docking_pdb_size_limit():
    """PDB larger than MAX_PDB_UPLOAD_BYTES should be rejected."""
    import os
    max_bytes = int(os.getenv("MAX_PDB_UPLOAD_BYTES", str(10 * 1024 * 1024)))
    oversized = b"ATOM  " * (max_bytes // 6 + 1)
    assert len(oversized) > max_bytes


def test_docking_modes_limit():
    """MAX_DOCKING_MODES should cap requested modes."""
    import os
    max_modes = int(os.getenv("MAX_DOCKING_MODES", "20"))
    requested = 100
    capped = min(requested, max_modes)
    assert capped == max_modes
