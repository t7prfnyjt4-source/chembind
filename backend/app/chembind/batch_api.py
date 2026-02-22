import os
import time
import uuid
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from app.chembind.celery_app import celery_app

router = APIRouter(prefix="/api", tags=["batch"])

# --- Firestore client (server-side; rules don't block service accounts) ---
_firestore_client = None

def _get_firestore():
    global _firestore_client
    if _firestore_client is not None:
        return _firestore_client
    from google.cloud import firestore  # type: ignore
    project_id = os.getenv("FIRESTORE_PROJECT_ID") or None
    _firestore_client = firestore.Client(project=project_id)
    return _firestore_client

def _jobs():
    return _get_firestore().collection("batch_jobs")

class BatchCreateRequest(BaseModel):
    rows: List[Dict[str, Any]] = Field(default_factory=list)
    attempt_cap: int = 3

class BatchCreateResponse(BaseModel):
    job_id: str

class BatchStatusResponse(BaseModel):
    job_id: str
    data: Dict[str, Any]

@router.post("/batch", response_model=BatchCreateResponse)
def create_batch(req: BatchCreateRequest) -> BatchCreateResponse:
    if not isinstance(req.rows, list) or len(req.rows) == 0:
        raise HTTPException(status_code=400, detail="rows must be a non-empty list")

    job_id = str(uuid.uuid4())

    # Create Firestore job doc (queued)
    _jobs().document(job_id).set(
        {
            "status": "queued",
            "attempt": 1,
            "attemptCap": int(req.attempt_cap),
            "total": len(req.rows),
            "processed": 0,
            "succeeded": 0,
            "failed": 0,
            "progress": 0.0,
            "createdAt": int(time.time()),
            "updatedAt": int(time.time()),
        },
        merge=True,
    )

    payload = {
        "job_id": job_id,
        "rows": req.rows,
        "attempt": 1,
        "attempt_cap": int(req.attempt_cap),
        "startedAt": int(time.time()),
    }

    celery_app.send_task("app.chembind.tasks.process_batch_job", args=[payload])

    return BatchCreateResponse(job_id=job_id)

@router.get("/batch/{job_id}", response_model=BatchStatusResponse)
def get_batch(job_id: str) -> BatchStatusResponse:
    doc = _jobs().document(job_id).get()
    if not doc.exists:
        raise HTTPException(status_code=404, detail="job not found")
    data = doc.to_dict() or {}
    return BatchStatusResponse(job_id=job_id, data=data)

@router.delete("/batch/{job_id}")
def delete_batch(job_id: str) -> Dict[str, Any]:
    # Optional cleanup endpoint (handy during testing)
    _jobs().document(job_id).delete()
    return {"ok": True, "job_id": job_id}
