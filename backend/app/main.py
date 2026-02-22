from __future__ import annotations

from pathlib import Path
from dotenv import load_dotenv

# Always load backend/.env
ENV_PATH = Path(__file__).resolve().parents[1] / ".env"
load_dotenv(dotenv_path=ENV_PATH)

import os
import uuid
from typing import Optional, Dict, Any

from app.chembind.jobs_api import router as jobs_router

from fastapi import FastAPI, Depends, Header, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from app.chembind.rdkit_safe import compute_descriptors, SmilesValidationError, RdkitLimits
from app.chembind.timeout_runner import run_with_timeout, TimeoutConfig, TimeoutError as HardTimeout
from app.chembind.firebase_admin import extract_bearer_token, verify_bearer_token
from app.chembind.firestore_repo import FirestoreRepo
from app.chembind.celery_app import celery_app  # adjust path if needed


# -------------------------
# Config (env driven)
# -------------------------
DEBUG = os.getenv("DEBUG", "0") == "1"
CORS_ORIGINS = [o.strip() for o in os.getenv("CORS_ORIGINS", "").split(",") if o.strip()]
ANALYZE_TIMEOUT_S = float(os.getenv("ANALYZE_TIMEOUT_S", "2.5"))
MAX_SMILES_LEN = int(os.getenv("MAX_SMILES_LEN", "500"))
MAX_ATOMS = int(os.getenv("MAX_ATOMS", "200"))

# Batch config
MAX_BATCH_ROWS = int(os.getenv("MAX_BATCH_ROWS", "500"))
REQUIRE_IDEMPOTENCY_KEY = os.getenv("REQUIRE_IDEMPOTENCY_KEY", "true").lower() == "true"
IDEMPOTENCY_TTL_HOURS = int(os.getenv("IDEMPOTENCY_TTL_HOURS", "24"))


# -------------------------
# Structured error response
# -------------------------
def err(code: str, message: str, status: int = 400):
    return JSONResponse(
        status_code=status,
        content={"error": {"code": code, "message": message}},
    )


# -------------------------
# Auth dependencies
# -------------------------
def get_optional_user(authorization: Optional[str] = Header(default=None)) -> Optional[Dict[str, Any]]:
    token = extract_bearer_token(authorization)
    if not token:
        return None
    try:
        decoded = verify_bearer_token(token)
        return {"uid": decoded.get("uid"), "claims": decoded}
    except Exception:
        return None


def get_required_user(user: Optional[Dict[str, Any]] = Depends(get_optional_user)) -> Dict[str, Any]:
    if not user or not user.get("uid"):
        raise HTTPException(status_code=401, detail="UNAUTHORIZED")
    return user


# -------------------------
# Schemas
# -------------------------
class AnalyzeRequest(BaseModel):
    smiles: str = Field(..., max_length=500)


class AnalyzeResponse(BaseModel):
    smiles: str
    descriptors: Dict[str, Any]


class BatchRow(BaseModel):
    smiles: str = Field(..., min_length=1)
    meta: Optional[Dict[str, Any]] = None


class CreateBatchRequest(BaseModel):
    source: str = Field("csv")
    rows: list[BatchRow]


class CreateBatchResponse(BaseModel):
    jobId: str
    status: str
    total: int
    idempotent: bool


# -------------------------
# App
# -------------------------
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

# -------------------------
# Routers
# -------------------------
from app.chembind.batch_api import router as batch_router
from app.chembind.jobs_api import router as jobs_router

app.include_router(batch_router)
app.include_router(jobs_router)

# -------------------------
# CORS
# -------------------------
if CORS_ORIGINS:
    app.add_middleware(
        CORSMiddleware,
        allow_origins=CORS_ORIGINS,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )
# -------------------------
# Exception handlers
# -------------------------
@app.exception_handler(HTTPException)
async def http_exception_handler(_, exc: HTTPException):
    if exc.status_code == 401:
        return err("UNAUTHORIZED", "Missing/invalid token", status=401)
    if isinstance(exc.detail, str) and exc.detail.isupper():
        return err(exc.detail, "Request failed", status=exc.status_code)
    return err("HTTP_ERROR", str(exc.detail), status=exc.status_code)


@app.exception_handler(SmilesValidationError)
async def smiles_handler(_, exc: SmilesValidationError):
    return err("INVALID_SMILES", str(exc), status=400)


@app.exception_handler(HardTimeout)
async def timeout_handler(_, exc: HardTimeout):
    return err("TIMEOUT", str(exc), status=408)


@app.exception_handler(Exception)
async def unhandled_handler(_, exc: Exception):
    if DEBUG:
        return err("INTERNAL_ERROR_DEBUG", repr(exc), status=500)
    return err("INTERNAL_ERROR", "Unexpected server error", status=500)


# -------------------------
# Routes
# -------------------------
@app.get("/api/health")
async def health():
    return {"ok": True}


@app.post("/api/analyze", response_model=AnalyzeResponse)
async def analyze(req: AnalyzeRequest, user: Optional[Dict[str, Any]] = Depends(get_optional_user)):
    limits = RdkitLimits(max_smiles_len=MAX_SMILES_LEN, max_atoms=MAX_ATOMS)
    cfg = TimeoutConfig(seconds=ANALYZE_TIMEOUT_S)

    descriptors = run_with_timeout(
        compute_descriptors,
        cfg,
        req.smiles,
        limits=limits,
    )

    if user:
        repo = FirestoreRepo()
        repo.save_analysis(
            user["uid"],
            {
                "smiles": req.smiles,
                "descriptors": descriptors,
            },
        )

    return {"smiles": req.smiles, "descriptors": descriptors}


@app.post("/api/batch", response_model=CreateBatchResponse)
async def create_batch(
    payload: CreateBatchRequest,
    user: Dict[str, Any] = Depends(get_required_user),
    idempotency_key: Optional[str] = Header(default=None, alias="Idempotency-Key"),
):
    uid = user["uid"]

    if REQUIRE_IDEMPOTENCY_KEY and not idempotency_key:
        raise HTTPException(status_code=400, detail="Idempotency-Key required")

    if payload.source != "csv":
        raise HTTPException(status_code=400, detail="source must be 'csv'")

    total = len(payload.rows)
    if total == 0:
        raise HTTPException(status_code=400, detail="rows cannot be empty")
    if total > MAX_BATCH_ROWS:
        raise HTTPException(status_code=400, detail=f"Too many rows (max {MAX_BATCH_ROWS})")

    repo = FirestoreRepo()

    # Idempotency check
    if idempotency_key:
        existing_job_id = repo.get_idempotent_job(uid, idempotency_key)
        if existing_job_id:
            return CreateBatchResponse(
                jobId=existing_job_id,
                status="queued",
                total=total,
                idempotent=True,
            )

    job_id = uuid.uuid4().hex

    repo.create_job(uid=uid, job_id=job_id, total=total, source=payload.source, max_attempts=3)

    for idx, row in enumerate(payload.rows):
        item_id = f"{idx:06d}"
        repo.write_job_item(
            uid=uid,
            job_id=job_id,
            item_id=item_id,
            payload={
                "index": idx,
                "input": row.dict(),
                "status": "pending",
            },
        )

    if idempotency_key:
        repo.set_idempotent_job(uid, idempotency_key, job_id, ttl_hours=IDEMPOTENCY_TTL_HOURS)

    celery_app.send_task("chembind.process_batch_job", args=[uid, job_id])

    return CreateBatchResponse(jobId=job_id, status="queued", total=total, idempotent=False)


@app.get("/api/me")
async def me(user: Dict[str, Any] = Depends(get_required_user)):
    return {"uid": user["uid"]}
