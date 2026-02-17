from typing import Optional, Dict, Any, List
from firebase_admin import firestore
from app.firebase_init import initialize_firebase

def _db():
    # Ensure Firebase Admin is initialized BEFORE creating firestore client
    initialize_firebase()
    return firestore.client()

def save_analysis(uid: str, analysis: Dict[str, Any]) -> str:
    db = _db()
    ref = db.collection("users").document(uid).collection("analyses").document()
    ref.set(analysis)
    return ref.id

def list_analyses(uid: str, limit: int = 20) -> List[Dict[str, Any]]:
    db = _db()
    docs = (
        db.collection("users")
        .document(uid)
        .collection("analyses")
        .order_by("createdAt", direction=firestore.Query.DESCENDING)
        .limit(limit)
        .stream()
    )
    out = []
    for d in docs:
        item = d.to_dict()
        item["id"] = d.id
        out.append(item)
    return out

def get_analysis(uid: str, analysis_id: str) -> Optional[Dict[str, Any]]:
    db = _db()
    doc = (
        db.collection("users")
        .document(uid)
        .collection("analyses")
        .document(analysis_id)
        .get()
    )
    if not doc.exists:
        return None
    item = doc.to_dict()
    item["id"] = doc.id
    return item



