from __future__ import annotations

from typing import Any, Dict, Optional
from fastapi import APIRouter, Depends, Header, HTTPException

from app.chembind.firebase_admin import extract_bearer_token, verify_bearer_token
from app.chembind.firestore_repo import FirestoreRepo

router = APIRouter(prefix="/api", tags=["jobs"])


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


@router.get("/jobs/{job_id}")
def get_job(job_id: str, user: Dict[str, Any] = Depends(get_required_user)):
    repo = FirestoreRepo()
    doc = repo.db.collection("users").document(user["uid"]).collection("jobs").document(job_id).get()
    if not doc.exists:
        raise HTTPException(status_code=404, detail="NOT_FOUND")
    out = doc.to_dict() or {}
    out["jobId"] = job_id
    return out


@router.get("/jobs/{job_id}/items")
def get_job_items(job_id: str, user: Dict[str, Any] = Depends(get_required_user), limit: int = 2000):
    repo = FirestoreRepo()
    items = repo.list_job_items(user["uid"], job_id, limit=limit)
    return {"items": items}
