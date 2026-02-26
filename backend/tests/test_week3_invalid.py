import pytest

INVALID_SMILES = [
    "INVALID",
    "XXYYZZ",
    "123456",
    "!@#$%^",
    "C(C(C(C(C(C(C",  # unclosed parentheses
    "C1CCCC",          # unclosed ring
    "Zebra",           # not SMILES
    "H2O",             # not SMILES
]

@pytest.mark.parametrize("s", INVALID_SMILES)
def test_invalid_smiles_returns_error_or_4xx(client, s):
    r = client.post("/api/analyze", json={"smiles": s})

    # Different implementations choose different contracts:
    # - 200 with {valid:false, error:...}
    # - 400/422 with {detail:...}
    assert r.status_code in (200, 400, 422)

    if r.status_code == 200:
        data = r.json()
        assert data.get("valid") is False
        # error may be string or structured
        assert ("error" in data) or ("detail" in data) or ("message" in data)


def test_500_char_garbage_rejected_or_invalid(client):
    s = "X" * 500
    r = client.post("/api/analyze", json={"smiles": s})
    assert r.status_code in (200, 400, 413, 422)

    if r.status_code == 200:
        data = r.json()
        assert data.get("valid") is False
