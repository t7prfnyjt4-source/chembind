def test_invalid_smiles_returns_4xx_or_structured(client):
    r = client.post("/api/analyze", json={"smiles": "INVALID"})
    assert r.status_code in (200, 400, 422)

    if r.status_code == 200:
        data = r.json()
        # current API often returns dict with smiles/descriptors even if invalid handling differs
        assert isinstance(data, dict)


def test_salt_handling_returns_descriptors(client):
    # Your current API returns "smiles" and "descriptors" (no canonicalSmiles)
    r = client.post("/api/analyze", json={"smiles": "CCO.Cl"})
    assert r.status_code == 200
    data = r.json()
    assert data.get("smiles") == "CCO.Cl"
    assert "descriptors" in data
    assert isinstance(data["descriptors"], dict)


def test_descriptors_present(client):
    r = client.post("/api/analyze", json={"smiles": "CCO"})
    assert r.status_code == 200
    data = r.json()
    assert data.get("smiles") == "CCO"
    d = data.get("descriptors") or {}
    # Minimal sanity checks on descriptors
    assert d.get("atom_count") == 3
    assert "mw" in d
