# app/chembind/firestore_repo.py
from __future__ import annotations

from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Optional

from google.cloud import firestore

from .firebase_admin import get_firestore_client


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _now_iso() -> str:
    return _now().isoformat()


class FirestoreRepo:
    def __init__(self):
        self.db = get_firestore_client()

    def list_job_items(self, uid: str, job_id: str, limit: int = 2000) -> List[Dict[str, Any]]:
        col = (
            self.db.collection("users")
            .document(uid)
            .collection("jobs")
            .document(job_id)
            .collection("items")
        )

        qs = col.order_by("created_at").limit(limit).stream()

        out: List[Dict[str, Any]] = []
        for d in qs:
            item = d.to_dict() or {}
            item["id"] = d.id
            item.setdefault("itemId", d.id)
            out.append(item)
        return out

    # -------------------------
    # Analyses
    # users/{uid}/analyses/{analysisId}
    # -------------------------
    def save_analysis(self, uid: str, record: Dict[str, Any]) -> str:
        col = self.db.collection("users").document(uid).collection("analyses")
        doc_ref = col.document()
        payload = dict(record)
        payload.setdefault("created_at", _now_iso())
        doc_ref.set(payload)
        return doc_ref.id

    def list_analyses(self, uid: str, limit: int = 50) -> List[Dict[str, Any]]:
        col = self.db.collection("users").document(uid).collection("analyses")
        qs = col.order_by("created_at", direction="DESCENDING").limit(limit).stream()
        out: List[Dict[str, Any]] = []
        for d in qs:
            item = d.to_dict() or {}
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
        item = doc.to_dict() or {}
        item["id"] = doc.id
        return item

    # -------------------------
    # Jobs
    # users/{uid}/jobs/{jobId}
    # users/{uid}/jobs/{jobId}/items/{itemId}
    # -------------------------
    def create_job(
        self,
        uid: str,
        job: Optional[Dict[str, Any]] = None,
        *,
        job_id: Optional[str] = None,
        total: Optional[int] = None,
        status: str = "queued",
        processed: int = 0,
        success_count: int = 0,
        failure_count: int = 0,
        attempt: int = 0,
        max_attempts: int = 3,
        source: Optional[str] = None,
    ) -> str:
        """
        Compatible with:
          - create_job(uid, job_dict) -> auto id
          - create_job(uid, job_id=..., total=..., ...) -> fixed id
        """
        col = self.db.collection("users").document(uid).collection("jobs")
        doc_ref = col.document(job_id) if job_id else col.document()

        payload: Dict[str, Any] = {}
        if job:
            payload.update(dict(job))

        # canonical fields
        if total is not None:
            payload["total"] = int(total)
        payload.setdefault("status", status)
        payload.setdefault("processed", int(processed))
        payload.setdefault("successCount", int(success_count))
        payload.setdefault("failureCount", int(failure_count))
        payload.setdefault("attempt", int(attempt))
        payload.setdefault("maxAttempts", int(max_attempts))
        if source is not None:
            payload["source"] = source

        payload.setdefault("created_at", _now_iso())
        payload["updated_at"] = _now_iso()
        payload["jobId"] = doc_ref.id

        doc_ref.set(payload, merge=True)
        return doc_ref.id

    def update_job(self, uid: str, job_id: str, patch: Dict[str, Any]) -> None:
        doc_ref = (
            self.db.collection("users")
            .document(uid)
            .collection("jobs")
            .document(job_id)
        )
        payload = dict(patch)
        payload["updated_at"] = _now_iso()
        doc_ref.set(payload, merge=True)

    def write_job_item(
        self,
        uid: str,
        job_id: str,
        item_id: str,
        payload: Optional[Dict[str, Any]] = None,
        data: Optional[Dict[str, Any]] = None,
        **_ignored: Any,
    ) -> None:
        """
        Accepts either payload=... or data=... (for compatibility).
        """
        if payload is None and data is not None:
            payload = data
        if payload is None:
            payload = {}

        doc_ref = (
            self.db.collection("users")
            .document(uid)
            .collection("jobs")
            .document(job_id)
            .collection("items")
            .document(item_id)
        )

        out = dict(payload)
        out.setdefault("created_at", _now_iso())
        out["updated_at"] = _now_iso()
        out.setdefault("itemId", item_id)

        doc_ref.set(out, merge=True)

    def list_job_items(self, uid: str, job_id: str, limit: int = 5000) -> List[Dict[str, Any]]:
        col = (
            self.db.collection("users")
            .document(uid)
            .collection("jobs")
            .document(job_id)
            .collection("items")
        )
        # If created_at is missing on some docs, order_by can fail.
        # We still try ordering; if it fails, fall back to unordered stream.
        try:
            qs = col.order_by("created_at", direction="ASCENDING").limit(limit).stream()
        except Exception:
            qs = col.limit(limit).stream()

        out: List[Dict[str, Any]] = []
        for d in qs:
            item = d.to_dict() or {}
            item["id"] = d.id
            out.append(item)
        return out

    # -------------------------
    # Idempotency
    # users/{uid}/idempotency/{key}
    # { jobId, createdAt, expiresAt }
    # -------------------------
    def get_idempotent_job(self, uid: str, key: str) -> Optional[str]:
        doc_ref = (
            self.db.collection("users")
            .document(uid)
            .collection("idempotency")
            .document(key)
        )
        snap = doc_ref.get()
        if not snap.exists:
            return None

        data = snap.to_dict() or {}
        job_id = data.get("jobId")
        expires_at = data.get("expiresAt")

        # expiresAt can be Firestore Timestamp or datetime
        try:
            exp_dt: Optional[datetime]
            if hasattr(expires_at, "to_datetime"):
                exp_dt = expires_at.to_datetime()
            elif isinstance(expires_at, datetime):
                exp_dt = expires_at
            else:
                exp_dt = None

            if exp_dt and exp_dt < _now():
                return None
        except Exception:
            pass

        return job_id if isinstance(job_id, str) and job_id else None

    def set_idempotent_job(self, uid: str, key: str, job_id: str, ttl_hours: int = 24) -> None:
        doc_ref = (
            self.db.collection("users")
            .document(uid)
            .collection("idempotency")
            .document(key)
        )
        expires_at = _now() + timedelta(hours=int(ttl_hours))

        doc_ref.set(
            {
                "jobId": job_id,
                "createdAt": firestore.SERVER_TIMESTAMP,
                "expiresAt": expires_at,
            },
            merge=True,
        )
