import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)

INVALID_SMILES = [
    "INVALID",
    "XXYYZZ",
    "123456",
    "!@#$%^",
    "C(C(C(C(C(C(C",  # unclosed parentheses
    "C1CCCC",          # unclosed ring
    "Zebra",
    "H2O",             # not SMILES
]

@pytest.mark.parametrize("s", INVALID_SMILES)
def test_invalid_smiles_structured_error(s):
    r = client.post("/analyze", json={"smiles": s})
    assert r.status_code == 200
    data = r.json()
    assert data["valid"] is False
    assert data["error"] is not None
    assert "Invalid SMILES" in data["error"]
    assert data["descriptors"] is None
    assert data["canonicalSmiles"] is None

def test_500_char_garbage():
    s = "X" * 500
    r = client.post("/analyze", json={"smiles": s})
    assert r.status_code == 200
    assert r.json()["valid"] is False
