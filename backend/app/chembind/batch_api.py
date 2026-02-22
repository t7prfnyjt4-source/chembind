# app/chembind/batch_api.py
from __future__ import annotations

import os
import uuid
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, Header, HTTPException
from pydantic import BaseModel, Field

from app.chembind.firebase_admin import extract_bearer_token, verify_bearer_token
from app.chembind.firestore_repo import FirestoreRepo
from app.chembind.celery_app import celery_app

router = APIRouter(prefix="/api", tags=["batch"])


# -------------------------
# Auth dependency (NO circular import with main.py)
# -------------------------
def get_required_user(authorization: Optional[str] = Header(default=None)) -> Dict[str, Any]:
    token = extract_bearer_token(authorization)
    if not token:
        raise HTTPException(status_code=401, detail="UNAUTHORIZED")
    try:
        decoded = verify_bearer_token(token)
        uid = decoded.get("uid") or decoded.get("user_id") or decoded.get("sub")
        if not uid:
            raise HTTPException(status_code=401, detail="UNAUTHORIZED")
        return {"uid": uid, "claims": decoded}
    except Exception:
        raise HTTPException(status_code=401, detail="UNAUTHORIZED")


# -------------------------
# Config
# -------------------------
MAX_BATCH_ROWS = int(os.getenv("MAX_BATCH_ROWS", "500"))
REQUIRE_IDEMPOTENCY_KEY = os.getenv("REQUIRE_IDEMPOTENCY_KEY", "true").lower() == "true"
IDEMPOTENCY_TTL_HOURS = int(os.getenv("IDEMPOTENCY_TTL_HOURS", "24"))


# -------------------------
# Schemas
# -------------------------
class BatchRow(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=500)


class BatchCreateRequest(BaseModel):
    source: str = Field(default="csv")
    rows: List[BatchRow] = Field(default_factory=list)


class BatchCreateResponse(BaseModel):
    jobId: str
    status: str
    total: int
    idempotent: bool


# -------------------------
# Routes
# -------------------------
@router.post("/batch", response_model=BatchCreateResponse)
def create_batch(
    req: BatchCreateRequest,
    user: Dict[str, Any] = Depends(get_required_user),
    idempotency_key: Optional[str] = Header(default=None, alias="Idempotency-Key"),
) -> BatchCreateResponse:
    uid = user["uid"]

    if not req.rows or len(req.rows) == 0:
        raise HTTPException(status_code=400, detail="rows must be a non-empty list")

    if len(req.rows) > MAX_BATCH_ROWS:
        raise HTTPException(status_code=400, detail=f"too many rows (max {MAX_BATCH_ROWS})")

    if REQUIRE_IDEMPOTENCY_KEY and not idempotency_key:
        raise HTTPException(status_code=400, detail="Idempotency-Key header required")

    repo = FirestoreRepo()

    # 1) Idempotency check
    if idempotency_key:
        existing_job_id = repo.get_idempotent_job(uid, idempotency_key)
        if existing_job_id:
            return BatchCreateResponse(
                jobId=existing_job_id,
                status="queued",
                total=len(req.rows),
                idempotent=True,
            )

    # 2) Create new job
    job_id = uuid.uuid4().hex
    total = len(req.rows)

    repo.create_job(
        uid=uid,
        job_id=job_id,
        total=total,
        status="queued",
        attempt=0,
        max_attempts=3,
        source=req.source,
    )

    # 3) Store input rows as items so the worker can load them
    # (tasks.py loads rows = repo.list_job_items(uid, job_id))
    for idx, row in enumerate(req.rows):
        repo.write_job_item(
            uid,
            job_id,
            item_id=str(idx),
            payload={
                "index": idx,
                "smiles": row.smiles.strip(),
                "kind": "input",
            },
        )

    # 4) Set idempotency mapping (TTL)
    if idempotency_key:
        repo.set_idempotent_job(uid, idempotency_key, job_id, ttl_hours=IDEMPOTENCY_TTL_HOURS)

    # 5) Enqueue task (task signature: process_batch_job(uid, job_id))
    celery_app.send_task("chembind.process_batch_job", args=[uid, job_id])

    return BatchCreateResponse(jobId=job_id, status="queued", total=total, idempotent=False)
