# app/chembind/firestore_repo.py
from __future__ import annotations

from typing import Any, Dict, List, Optional
from datetime import datetime, timezone

from .firebase_admin import get_firestore_client


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


class FirestoreRepo:
    def __init__(self):
        self.db = get_firestore_client()

    # -------------------------
    # Analyses
    # users/{uid}/analyses/{analysisId}
    # -------------------------
    def save_analysis(self, uid: str, record: Dict[str, Any]) -> str:
        col = self.db.collection("users").document(uid).collection("analyses")
        doc_ref = col.document()
        record = dict(record)
        record.setdefault("created_at", _now_iso())
        doc_ref.set(record)
        return doc_ref.id

    def list_analyses(self, uid: str, limit: int = 50) -> List[Dict[str, Any]]:
        col = self.db.collection("users").document(uid).collection("analyses")
        qs = col.order_by("created_at", direction="DESCENDING").limit(limit).stream()
        out = []
        for d in qs:
            item = d.to_dict()
            item["id"] = d.id
            out.append(item)
        return out

    def get_analysis(self, uid: str, analysis_id: str) -> Optional[Dict[str, Any]]:
        doc = (
            self.db.collection("users")
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

    # -------------------------
    # Jobs (placeholder for Segment 2/3)
    # users/{uid}/jobs/{jobId}
    # users/{uid}/jobs/{jobId}/items/{itemId}
    # users/{uid}/idempotency/{key}
    # -------------------------
    def create_job(self, uid: str, job: Dict[str, Any]) -> str:
        col = self.db.collection("users").document(uid).collection("jobs")
        doc_ref = col.document()
        job = dict(job)
        job.setdefault("created_at", _now_iso())
        job.setdefault("status", "queued")
        doc_ref.set(job)
        return doc_ref.id

    def add_job_item(self, uid: str, job_id: str, item: Dict[str, Any]) -> str:
        col = (
            self.db.collection("users")
            .document(uid)
            .collection("jobs")
            .document(job_id)
            .collection("items")
        )
        doc_ref = col.document()
        item = dict(item)
        item.setdefault("created_at", _now_iso())
        doc_ref.set(item)
        return doc_ref.id

    def set_idempotency_key(self, uid: str, key: str, payload: Dict[str, Any]) -> None:
        doc = self.db.collection("users").document(uid).collection("idempotency").document(key)
        payload = dict(payload)
        payload.setdefault("created_at", _now_iso())
        doc.set(payload)

    def get_idempotency_key(self, uid: str, key: str) -> Optional[Dict[str, Any]]:
        doc = self.db.collection("users").document(uid).collection("idempotency").document(key).get()
        if not doc.exists:
            return None
        out = doc.to_dict()
        out["id"] = doc.id
        return out
