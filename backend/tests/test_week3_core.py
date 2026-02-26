import pytest

REFERENCE_MOLECULES = [
    {
        "name": "Ethanol",
        "smiles": "CCO",
        "mw": 46.07,
        "logp": -0.03,
        "hbd": 1,
        "hba": 1,
        "tpsa": 20.23,
        "atom_count": 3,
    },
    {
        "name": "Benzene",
        "smiles": "c1ccccc1",
        "mw": 78.11,
        "logp": 1.69,
        "hbd": 0,
        "hba": 0,
        "tpsa": 0.0,
        "atom_count": 6,
    },
    {
        "name": "Acetic Acid",
        "smiles": "CC(=O)O",
        "mw": 60.05,
        "logp": 0.0909,
        "hbd": 1,
        "hba": 1,
        "tpsa": 37.30,
        "atom_count": 4,
    },
    {
        "name": "Aspirin",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "mw": 180.16,
        "logp": 1.19,
        "hbd": 1,
        "hba": 3,
        "tpsa": 63.60,
        "atom_count": 13,
    },
]


@pytest.mark.parametrize("m", REFERENCE_MOLECULES)
def test_reference_molecule_descriptors(client, m):
    r = client.post("/api/analyze", json={"smiles": m["smiles"]})
    assert r.status_code == 200
    data = r.json()

    assert data.get("smiles") == m["smiles"]
    d = data.get("descriptors") or {}
    assert isinstance(d, dict)

    # numeric checks (tolerant)
    assert abs(d["mw"] - m["mw"]) < m["mw"] * 0.02
    assert abs(d["logp"] - m["logp"]) < 0.5
    assert d["hbd"] == m["hbd"]
    assert d["hba"] == m["hba"]
    assert d["atom_count"] == m["atom_count"]
    assert abs(d["tpsa"] - m["tpsa"]) < 2.0


def test_single_atom_methane(client):
    r = client.post("/api/analyze", json={"smiles": "C"})
    assert r.status_code == 200
    d = r.json()["descriptors"]
    assert d["atom_count"] == 1


def test_aromatic_notation(client):
    r = client.post("/api/analyze", json={"smiles": "c1ccccc1"})
    assert r.status_code == 200
    d = r.json()["descriptors"]
    assert d["atom_count"] == 6


def test_charged_species_ammonium(client):
    r = client.post("/api/analyze", json={"smiles": "[NH4+]"})
    assert r.status_code == 200
    d = r.json()["descriptors"]
    assert d["atom_count"] == 1


def test_empty_string_rejected(client):
    r = client.post("/api/analyze", json={"smiles": ""})
    assert r.status_code in (400, 422)


def test_whitespace_only_rejected(client):
    r = client.post("/api/analyze", json={"smiles": "   "})
    assert r.status_code in (400, 422)


def test_salt_multiple_fragments_still_returns_descriptors(client):
    r = client.post("/api/analyze", json={"smiles": "CCO.Cl.[Na+]"} )
    assert r.status_code == 200
    data = r.json()
    assert data.get("smiles") == "CCO.Cl.[Na+]"
    assert "descriptors" in data
