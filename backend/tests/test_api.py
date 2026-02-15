from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)

def test_root_ok():
    r = client.get("/")
    assert r.status_code == 200
    assert r.json().get("ok") is True

def test_analyze_ethanol():
    r = client.post("/analyze", json={"smiles": "CCO"})
    assert r.status_code == 200
    data = r.json()
    assert data["valid"] is True
    assert data["canonicalSmiles"] == "CCO"
    for k in ["mw", "logp", "hbd", "hba", "tpsa", "atom_count"]:
        assert k in data["descriptors"]
