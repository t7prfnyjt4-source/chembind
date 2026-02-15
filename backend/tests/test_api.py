from fastapi.testclient import TestClient
from main import app

client = TestClient(app)

def test_root_ok():
    r = client.get("/")
    assert r.status_code == 200
    assert r.json().get("ok") is True

def test_analyze_ethanol():
    r = client.get("/analyze", params={"smiles": "CCO"})
    assert r.status_code == 200
    data = r.json()
    assert data["smiles"] == "CCO"
    for k in ["mw", "logp", "hbd", "hba", "tpsa", "num_atoms"]:
        assert k in data
