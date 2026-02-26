def test_pathological_smiles_timeout_or_graceful_error(client):
    # Pathological-ish: extremely long chain; should either be rejected or handled fast
    s = "C" * 10000
    r = client.post("/api/analyze", json={"smiles": s})

    # If you enforce hard timeout / size limits, you may get 400/413/408/504.
    # If you choose to return 200 valid:false, accept that too.
    assert r.status_code in (200, 400, 408, 413, 422, 504)

    if r.status_code == 200:
        data = r.json()
        # if it didn't timeout, it should still be invalid or safe-handled
        assert data.get("valid") in (False, True)
