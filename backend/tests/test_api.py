def test_root_ok(client):
    r = client.get("/api/health")
    assert r.status_code == 200
    assert r.json().get("ok") is True


def test_analyze_ethanol(client):
    r = client.post("/api/analyze", json={"smiles": "CCO"})
    assert r.status_code == 200
    data = r.json()
    # keep assertions minimal; downstream schema may evolve
    assert "canonicalSmiles" in data or "result" in data or "descriptors" in data
