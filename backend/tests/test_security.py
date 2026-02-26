def test_security_headers_present_on_health(client):
    r = client.get("/api/health")
    assert r.status_code == 200

    # headers you confirmed exist in backend/app/main.py
    assert r.headers.get("x-content-type-options") == "nosniff"
    assert r.headers.get("x-frame-options") == "DENY"
    assert r.headers.get("referrer-policy") == "no-referrer"
    assert "geolocation" in (r.headers.get("permissions-policy") or "")


def test_security_headers_present_on_analyze(client):
    r = client.post("/api/analyze", json={"smiles": "CCO"})
    assert r.status_code == 200

    assert r.headers.get("x-content-type-options") == "nosniff"
    assert r.headers.get("x-frame-options") == "DENY"
    assert r.headers.get("referrer-policy") == "no-referrer"
    assert "geolocation" in (r.headers.get("permissions-policy") or "")
