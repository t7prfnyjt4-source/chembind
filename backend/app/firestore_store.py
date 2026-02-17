from __future__ import annotations

from typing import Any, Dict, List, Optional
from firebase_admin import firestore

def _db():
    # IMPORTANT: firestore.client() requires firebase_admin.initialize_app() to have already run.
    return firestore.client()

def save_analysis(uid: str, doc: Dict[str, Any]) -> str:
    ref = _db().collection("users").document(uid).collection("analyses").document()
    ref.set(doc)
    return ref.id

def list_analyses(uid: str, limit: int = 50) -> List[Dict[str, Any]]:
    q = (
        _db()
        .collection("users").document(uid)
        .collection("analyses")
        .order_by("createdAt", direction=firestore.Query.DESCENDING)
        .limit(limit)
    )
    return [{"id": d.id, **d.to_dict()} for d in q.stream()]

def get_analysis(uid: str, analysis_id: str) -> Optional[Dict[str, Any]]:
    snap = (
        _db()
        .collection("users").document(uid)
        .collection("analyses")
        .document(analysis_id)
        .get()
    )
    if not snap.exists:
        return None
    return {"id": snap.id, **snap.to_dict()}




