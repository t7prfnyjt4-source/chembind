import os
import uuid
import pytest

def _auth_headers():
    """
    Provide a token for local tests via env var:
      export CHEMBIND_TEST_BEARER="YOUR_TOKEN"
    """
    token = os.getenv("CHEMBIND_TEST_BEARER", "").strip()
    if not token:
        return None
    return {"Authorization": f"Bearer {token}"}


def test_batch_requires_auth(client):
    r = client.post("/api/batch", json={"rows": [{"smiles": "CCO"}]})
    assert r.status_code in (401, 403)


@pytest.mark.skipif(_auth_headers() is None, reason="Set CHEMBIND_TEST_BEARER to run authenticated batch tests")
def test_batch_valid_minimal_authed(client):
    headers = _auth_headers()
    r = client.post("/api/batch", json={"rows": [{"smiles": "CCO"}]}, headers=headers)
    assert r.status_code in (200, 201, 202), r.text[:200]
    data = r.json()
    assert any(k in data for k in ("jobId", "job_id", "id")), data


@pytest.mark.skipif(_auth_headers() is None, reason="Set CHEMBIND_TEST_BEARER to run authenticated batch tests")
def test_batch_empty_rows_rejected_authed(client):
    headers = _auth_headers()
    r = client.post("/api/batch", json={"rows": []}, headers=headers)
    assert r.status_code in (400, 422), r.text[:200]


@pytest.mark.skipif(_auth_headers() is None, reason="Set CHEMBIND_TEST_BEARER to run authenticated batch tests")
def test_batch_idempotency_key_dedup_authed(client):
    headers = _auth_headers()
    key = str(uuid.uuid4())

    h1 = dict(headers)
    h1["Idempotency-Key"] = key

    r1 = client.post("/api/batch", json={"rows": [{"smiles": "CCO"}]}, headers=h1)
    r2 = client.post("/api/batch", json={"rows": [{"smiles": "CCO"}]}, headers=h1)

    assert r1.status_code in (200, 201, 202), r1.text[:200]
    assert r2.status_code in (200, 201, 202), r2.text[:200]

    d1 = r1.json()
    d2 = r2.json()
    id1 = d1.get("jobId") or d1.get("job_id") or d1.get("id")
    id2 = d2.get("jobId") or d2.get("job_id") or d2.get("id")
    assert id1 == id2
