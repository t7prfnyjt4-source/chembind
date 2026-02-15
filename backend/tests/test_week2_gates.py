from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)

def test_invalid_smiles():
    r = client.post("/analyze", json={"smiles": "INVALID"})
    data = r.json()
    assert data["valid"] is False
    assert "error" in data

def test_salt_handling():
    r = client.post("/analyze", json={"smiles": "CCO.Cl"})
    assert r.json()["canonicalSmiles"] == "CCO"

def test_canonical():
    r1 = client.post("/analyze", json={"smiles": "CCO"})
    r2 = client.post("/analyze", json={"smiles": "C(O)C"})
    assert r1.json()["canonicalSmiles"] == r2.json()["canonicalSmiles"]

def test_rounding():
    r = client.post("/analyze", json={"smiles": "CCO"})
    mw = str(r.json()["descriptors"]["mw"])
    assert len(mw.split(".")[-1]) <= 4
