import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)

REFERENCE_MOLECULES = [
    # values are approximate; tolerances are used
    {"name": "Ethanol", "smiles": "CCO", "canonical": "CCO", "mw": 46.07, "logp": -0.03, "hbd": 1, "hba": 1, "tpsa": 20.23, "atom_count": 3},
    {"name": "Benzene", "smiles": "c1ccccc1", "canonical": "c1ccccc1", "mw": 78.11, "logp": 1.69, "hbd": 0, "hba": 0, "tpsa": 0.0, "atom_count": 6},
    {"name": "Acetic Acid", "smiles": "CC(=O)O", "canonical": "CC(=O)O", "mw": 60.05, "logp": 0.0909, "hbd": 1, "hba": 1, "tpsa": 37.30, "atom_count": 4},
    {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "canonical": "CC(=O)Oc1ccccc1C(=O)O", "mw": 180.16, "logp": 1.19, "hbd": 1, "hba": 3, "tpsa": 63.60, "atom_count": 13},
]

@pytest.mark.parametrize("m", REFERENCE_MOLECULES)
def test_reference_molecule(m):
    r = client.post("/analyze", json={"smiles": m["smiles"]})
    assert r.status_code == 200
    data = r.json()
    assert data["valid"] is True
    assert data["canonicalSmiles"] == m["canonical"]

    d = data["descriptors"]
    assert abs(d["mw"] - m["mw"]) < m["mw"] * 0.01
    assert abs(d["logp"] - m["logp"]) < 0.35
    assert d["hbd"] == m["hbd"]
    assert d["hba"] == m["hba"]
    assert d["atom_count"] == m["atom_count"]
    assert abs(d["tpsa"] - m["tpsa"]) < 1.0

def test_single_atom_methane():
    r = client.post("/analyze", json={"smiles": "C"})
    data = r.json()
    assert data["valid"] is True
    assert data["descriptors"]["atom_count"] == 1

def test_aromatic_notation():
    r = client.post("/analyze", json={"smiles": "c1ccccc1"})
    data = r.json()
    assert data["valid"] is True
    assert data["canonicalSmiles"] == "c1ccccc1"

def test_charged_species_ammonium():
    r = client.post("/analyze", json={"smiles": "[NH4+]"})
    data = r.json()
    assert data["valid"] is True
    assert data["descriptors"]["atom_count"] == 1

def test_empty_string_rejected_422():
    r = client.post("/analyze", json={"smiles": ""})
    assert r.status_code == 422

def test_whitespace_only_rejected_422():
    r = client.post("/analyze", json={"smiles": "   "})
    assert r.status_code == 422

def test_salt_multiple_fragments_keeps_largest():
    r = client.post("/analyze", json={"smiles": "CCO.Cl.[Na+]"})
    data = r.json()
    assert data["valid"] is True
    assert data["canonicalSmiles"] == "CCO"

def test_canonical_variants_match():
    variants = ["CCO", "C(O)C", "OCC"]
    canon = []
    for s in variants:
        r = client.post("/analyze", json={"smiles": s})
        data = r.json()
        assert data["valid"] is True
        canon.append(data["canonicalSmiles"])
    assert len(set(canon)) == 1

def test_float_precision_max_4_decimals():
    r = client.post("/analyze", json={"smiles": "CCO"})
    d = r.json()["descriptors"]
    for k in ["mw", "logp", "tpsa"]:
        s = str(d[k])
        if "." in s:
            assert len(s.split(".")[-1]) <= 4
