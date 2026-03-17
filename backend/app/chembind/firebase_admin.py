# app/chembind/firebase_admin.py
from __future__ import annotations

import json
import os
from typing import Optional, Dict, Any

import firebase_admin
from firebase_admin import credentials, auth, firestore


def init_firebase() -> None:
    """
    Initialize Firebase Admin once.
    Uses:
    - FIREBASE_SERVICE_ACCOUNT_JSON (raw JSON string) — preferred for Railway/Render
    - GOOGLE_APPLICATION_CREDENTIALS (ADC) as fallback
    """
    if firebase_admin._apps:
        return

    sa_json = os.getenv("FIREBASE_SERVICE_ACCOUNT_JSON")
    if sa_json:
        cred = credentials.Certificate(json.loads(sa_json))
        firebase_admin.initialize_app(cred)
        return

    # fallback: ADC
    firebase_admin.initialize_app()


def get_firestore_client():
    init_firebase()
    return firestore.client()


def verify_bearer_token(id_token: str) -> Dict[str, Any]:
    init_firebase()
    decoded = auth.verify_id_token(id_token)
    return decoded


def extract_bearer_token(auth_header: Optional[str]) -> Optional[str]:
    if not auth_header:
        return None
    parts = auth_header.split()
    if len(parts) == 2 and parts[0].lower() == "bearer":
        return parts[1].strip()
    return None
