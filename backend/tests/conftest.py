import pytest
from fastapi.testclient import TestClient
from app.main import app

@pytest.fixture(scope="session")
def client():
    # Fail-safe: avoids TrustedHost 400 by making Host=localhost instead of testserver
    return TestClient(app, base_url="http://localhost")
