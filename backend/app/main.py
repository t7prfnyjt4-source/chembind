from __future__ import annotations

import os
from typing import Optional, Dict, Any

from fastapi import FastAPI, Depends, Header, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from app.chembind.rdkit_safe import compute_descriptors, SmilesValidationError, RdkitLimits
from app.chembind.timeout_runner import run_with_timeout, TimeoutConfig, TimeoutError as HardTimeout
from app.chembind.firebase_admin import extract_bearer_token, verify_bearer_token
from app.chembind.firestore_repo import FirestoreRepo


# -------------------------
# Config (env driven)
# -------------------------
DEBUG = os.getenv("DEBUG", "0") == "1"
CORS_ORIGINS = [o.strip() for o in os.getenv("CORS_ORIGINS", "").split(",") if o.strip()]
ANALYZE_TIMEOUT_S = float(os.getenv("ANALYZE_TIMEOUT_S", "2.5"))
MAX_SMILES_LEN = int(os.getenv("MAX_SMILES_LEN", "500"))
MAX_ATOMS = int(os.getenv("MAX_ATOMS", "200"))


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


# -------------------------
# App
# -------------------------
app = FastAPI()

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
# Routes (ALL under /api)
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


@app.get("/api/me")
async def me(user: Dict[str, Any] = Depends(get_required_user)):
    return {"uid": user["uid"]}


@app.get("/api/analyses")
async def list_analyses(user: Dict[str, Any] = Depends(get_required_user), limit: int = 50):
    repo = FirestoreRepo()
    return {"items": repo.list_analyses(user["uid"], limit=limit)}


@app.get("/api/analyses/{analysis_id}")
async def get_analysis(analysis_id: str, user: Dict[str, Any] = Depends(get_required_user)):
    repo = FirestoreRepo()
    item = repo.get_analysis(user["uid"], analysis_id)
    if not item:
        return err("NOT_FOUND", "Analysis not found", status=404)
    return item
