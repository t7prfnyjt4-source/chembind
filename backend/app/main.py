# backend/app/main.py

from __future__ import annotations

import logging
import os
import uuid
from typing import Optional, Dict, Any

from fastapi import FastAPI, Depends, Header, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.middleware.trustedhost import TrustedHostMiddleware

from app.chembind.ratelimit import build_limiter
from app.chembind.logger import request_id_var, setup_json_logging
from app.chembind.rdkit_safe import compute_descriptors, compute_morgan_fp, smiles_to_mol, SmilesValidationError, RdkitLimits
from app.chembind.similarity import tanimoto_search, substructure_search
from app.chembind.conformers import generate_conformers, conformer_cache_key
from app.chembind.timeout_runner import run_with_timeout, TimeoutConfig, TimeoutError as HardTimeout
from app.chembind.firebase_admin import extract_bearer_token, verify_bearer_token
from app.chembind.firestore_repo import FirestoreRepo
from app.chembind.celery_app import celery_app

from app.chembind.batch_api import router as batch_router
from app.chembind.jobs_api import router as jobs_router


# -------------------------
# Config (env driven)
# -------------------------
DEBUG = os.getenv("DEBUG", "0") == "1"

CORS_ORIGINS = [o.strip() for o in os.getenv("CORS_ORIGINS", "").split(",") if o.strip()]
ANALYZE_TIMEOUT_S = float(os.getenv("ANALYZE_TIMEOUT_S", "2.5"))

TRUSTED_HOSTS = [
    h.strip()
    for h in os.getenv("TRUSTED_HOSTS", "localhost,127.0.0.1").split(",")
    if h.strip()
]

MAX_REQUEST_BYTES = int(os.getenv("MAX_REQUEST_BYTES", str(2 * 1024 * 1024)))  # 2MB default
MAX_SMILES_LEN = int(os.getenv("MAX_SMILES_LEN", "500"))
MAX_ATOMS = int(os.getenv("MAX_ATOMS", "200"))

# Batch config
MAX_BATCH_ROWS = int(os.getenv("MAX_BATCH_ROWS", "500"))
REQUIRE_IDEMPOTENCY_KEY = os.getenv("REQUIRE_IDEMPOTENCY_KEY", "true").lower() == "true"
IDEMPOTENCY_TTL_HOURS = int(os.getenv("IDEMPOTENCY_TTL_HOURS", "24"))

# Feature flags
ENABLE_SIMILARITY_SEARCH = os.getenv("ENABLE_SIMILARITY_SEARCH", "false").lower() == "true"
ENABLE_CONFORMERS = os.getenv("ENABLE_CONFORMERS", "false").lower() == "true"
ENABLE_DOCKING = os.getenv("ENABLE_DOCKING", "false").lower() == "true"
ENABLE_EXPORT = os.getenv("ENABLE_EXPORT", "false").lower() == "true"
ENABLE_ANNOTATIONS = os.getenv("ENABLE_ANNOTATIONS", "false").lower() == "true"

# -------------------------
# Segment 6 — Structured logging + Sentry (env driven)
# -------------------------
setup_json_logging(os.getenv("LOG_LEVEL", "INFO"))
logger = logging.getLogger("chembind")

SENTRY_DSN = os.getenv("SENTRY_DSN", "").strip()
if SENTRY_DSN:
    import sentry_sdk
    sentry_sdk.init(
        dsn=SENTRY_DSN,
        environment=os.getenv("ENVIRONMENT", "dev"),
        traces_sample_rate=float(os.getenv("SENTRY_TRACES_SAMPLE_RATE", "0.0")),
    )
    logger.info("Sentry enabled")
else:
    logger.info("Sentry disabled (no SENTRY_DSN)")


# -------------------------
# Structured error response
# -------------------------
def err(code: str, message: str, status: int = 400):
    return JSONResponse(
        status_code=status,
        content={"error": {"code": code, "message": message}},
    )


# -------------------------
# Segment 4 — Security middlewares
# -------------------------
class RequestIdMiddleware(BaseHTTPMiddleware):
    """
    - Reads incoming X-Request-Id if present, else generates one.
    - Stores it in ContextVar for logging.
    - Adds X-Request-Id to every response.
    """

    async def dispatch(self, request, call_next):
        incoming = request.headers.get("x-request-id")
        rid = incoming.strip() if incoming else str(uuid.uuid4())

        token = request_id_var.set(rid)
        try:
            response = await call_next(request)
        finally:
            request_id_var.reset(token)

        response.headers["X-Request-Id"] = rid
        return response


class RequestSizeLimitMiddleware(BaseHTTPMiddleware):
    """
    Enforces MAX_REQUEST_BYTES against:
    - Content-Length header (if present)
    - Actual body size (always checked)
    """

    def __init__(self, app, max_bytes: int):
        super().__init__(app)
        self.max_bytes = max_bytes

    async def dispatch(self, request, call_next):
        cl = request.headers.get("content-length")
        if cl is not None:
            try:
                if int(cl) > self.max_bytes:
                    return JSONResponse(status_code=413, content={"detail": "Request too large"})
            except ValueError:
                return JSONResponse(status_code=400, content={"detail": "Invalid Content-Length header"})

        body = await request.body()
        if len(body) > self.max_bytes:
            return JSONResponse(status_code=413, content={"detail": "Request too large"})

        # Re-inject body for downstream handlers
        async def receive():
            return {"type": "http.request", "body": body, "more_body": False}

        request._receive = receive  # Starlette internal; common pattern
        return await call_next(request)


class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """
    OWASP-ish baseline response headers.
    """

    async def dispatch(self, request, call_next):
        response = await call_next(request)

        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["Permissions-Policy"] = "geolocation=(), microphone=(), camera=()"

        # HSTS only when HTTPS
        if request.url.scheme == "https":
            response.headers["Strict-Transport-Security"] = "max-age=31536000; includeSubDomains"

        return response


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
app = FastAPI()

# -----------------------------
# Segment 5 — Rate Limiter
# -----------------------------
log = logging.getLogger(__name__)
limiter = build_limiter()

# SlowAPI middleware expects a real slowapi.Limiter in app.state.limiter.
# Our decorators use `limiter` (SafeLimiter/NoopLimiter), but middleware must see the inner limiter when present.
try:
    app.state.limiter = limiter.inner()  # SafeLimiter -> real limiter
except AttributeError:
    app.state.limiter = limiter          # NoopLimiter or already-real limiter


# -----------------------------
# SlowAPI Middleware + 429
# -----------------------------
try:
    from slowapi.errors import RateLimitExceeded
    from slowapi.middleware import SlowAPIMiddleware
    from slowapi import _rate_limit_exceeded_handler

    # Only enable middleware when we have a real limiter (SafeLimiter wraps one)
    if hasattr(limiter, "inner"):
        app.add_middleware(SlowAPIMiddleware)
        app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

except Exception as e:
    log.warning(f"SlowAPI not enabled: {e}")

# ---------------------------
# Segment 4 middleware registration
# ---------------------------
app.add_middleware(TrustedHostMiddleware, allowed_hosts=TRUSTED_HOSTS)
app.add_middleware(RequestIdMiddleware)
app.add_middleware(RequestSizeLimitMiddleware, max_bytes=MAX_REQUEST_BYTES)
app.add_middleware(SecurityHeadersMiddleware)

# CORS
if CORS_ORIGINS:
    app.add_middleware(
        CORSMiddleware,
        allow_origins=CORS_ORIGINS,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

# Routers
# app.include_router(batch_router)
app.include_router(jobs_router)

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
    # Don’t leak internals unless DEBUG is explicitly enabled
    if DEBUG:
        return err("INTERNAL_ERROR_DEBUG", repr(exc), status=500)
    return err("INTERNAL_ERROR", "Unexpected server error", status=500)


# -------------------------
# Routes
# -------------------------
@app.get("/api/health")
async def health():
    return {"ok": True}


@app.head("/api/health")
async def health_head():
    return {}


# ✅ Replace BOTH @app.post(...) lines with this one

@app.post(
    "/api/analyze",
    response_model=AnalyzeResponse,
    dependencies=[Depends(limiter)],
)
async def analyze(
    request: Request,
    req: AnalyzeRequest,
    user: Optional[Dict[str, Any]] = Depends(get_optional_user),
):
    limits = RdkitLimits(max_smiles_len=MAX_SMILES_LEN, max_atoms=MAX_ATOMS)
    cfg = TimeoutConfig(seconds=ANALYZE_TIMEOUT_S)

    descriptors = run_with_timeout(
        compute_descriptors,
        cfg,
        req.smiles,
        limits=limits,
    )

    # Compute Morgan fingerprint (non-breaking: runs outside timeout for simplicity)
    morgan_fp: str | None = None
    try:
        mol = smiles_to_mol(req.smiles, limits)
        morgan_fp = compute_morgan_fp(mol)
    except Exception:
        pass  # fingerprint is optional; don't fail the request

    if user:
        repo = FirestoreRepo()
        record: dict = {"smiles": req.smiles, "descriptors": descriptors}
        if morgan_fp:
            record["morganFp2048"] = morgan_fp
        repo.save_analysis(user["uid"], record)

    return {"smiles": req.smiles, "descriptors": descriptors}


@app.post(
    "/api/batch",
    response_model=CreateBatchResponse,
    dependencies=[Depends(limiter)],
)
async def create_batch(
    request: Request,
    payload: CreateBatchRequest,
    user: Optional[Dict[str, Any]] = Depends(get_optional_user),
    idempotency_key: Optional[str] = Header(default=None, alias="Idempotency-Key"),
):

    # Enforce auth AFTER limiter runs
    if not user or not user.get("uid"):
        raise HTTPException(status_code=401, detail="UNAUTHORIZED")

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


# -------------------------
# Segment 23 — Similarity Search
# -------------------------
class SimilaritySearchRequest(BaseModel):
    smiles: str = Field(..., max_length=500)
    top_k: int = Field(10, ge=1, le=50)
    min_sim: float = Field(0.3, ge=0.0, le=1.0)


@app.post(
    "/api/similarity/search",
    dependencies=[Depends(limiter)],
)
async def similarity_search(
    request: Request,
    req: SimilaritySearchRequest,
    user: Dict[str, Any] = Depends(get_required_user),
):
    if not ENABLE_SIMILARITY_SEARCH:
        raise HTTPException(status_code=404, detail="Feature not enabled")

    repo = FirestoreRepo()
    history = repo.list_analyses(user["uid"], limit=500)

    results = tanimoto_search(
        query_smiles=req.smiles,
        history_docs=history,
        top_k=req.top_k,
        min_sim=req.min_sim,
    )

    return {"results": results}


# -------------------------
# Segment 24 — Unified Search (similarity + substructure)
# -------------------------
class UnifiedSearchRequest(BaseModel):
    query: str = Field(..., max_length=500)
    mode: str = Field("similarity")  # "similarity" or "substructure"
    top_k: int = Field(10, ge=1, le=50)
    min_sim: float = Field(0.3, ge=0.0, le=1.0)
    max_results: int = Field(50, ge=1, le=200)


@app.post(
    "/api/search",
    dependencies=[Depends(limiter)],
)
async def unified_search(
    request: Request,
    req: UnifiedSearchRequest,
    user: Dict[str, Any] = Depends(get_required_user),
):
    if not ENABLE_SIMILARITY_SEARCH:
        raise HTTPException(status_code=404, detail="Feature not enabled")

    if req.mode not in ("similarity", "substructure"):
        raise HTTPException(status_code=400, detail="mode must be 'similarity' or 'substructure'")

    repo = FirestoreRepo()
    history = repo.list_analyses(user["uid"], limit=500)

    if req.mode == "similarity":
        results = tanimoto_search(
            query_smiles=req.query,
            history_docs=history,
            top_k=req.top_k,
            min_sim=req.min_sim,
        )
    else:
        results = substructure_search(
            smarts=req.query,
            history_docs=history,
            max_results=req.max_results,
        )

    return {"results": results, "mode": req.mode}


# -------------------------
# Segment 26 — Conformer Ensemble
# -------------------------
def _get_redis():
    """Return a Redis client or None if unavailable."""
    try:
        import redis
        url = os.getenv("REDIS_URL", "redis://localhost:6379/0")
        return redis.from_url(url, decode_responses=True)
    except Exception:
        return None


class ConformerRequest(BaseModel):
    smiles: str = Field(..., max_length=500)
    num_confs: int = Field(20, ge=1, le=50)


@app.post(
    "/api/conformers",
    dependencies=[Depends(limiter)],
)
async def conformers_endpoint(
    request: Request,
    req: ConformerRequest,
    user: Dict[str, Any] = Depends(get_required_user),
):
    if not ENABLE_CONFORMERS:
        raise HTTPException(status_code=404, detail="Feature not enabled")

    import json as _json

    cache_key = conformer_cache_key(req.smiles, req.num_confs)
    rds = _get_redis()

    # Check cache
    if rds:
        try:
            cached = rds.get(cache_key)
            if cached:
                return {"conformers": _json.loads(cached), "cached": True}
        except Exception:
            pass

    confs = generate_conformers(
        smiles=req.smiles,
        num_confs=req.num_confs,
    )

    # Store in cache (1h TTL)
    if rds:
        try:
            rds.setex(cache_key, 3600, _json.dumps(confs))
        except Exception:
            pass

    return {"conformers": confs, "cached": False}
