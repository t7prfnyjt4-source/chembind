import os
import time
from typing import Any, Dict, List

from celery.exceptions import SoftTimeLimitExceeded
from celery.utils.log import get_task_logger

from app.chembind.celery_app import celery_app

log = get_task_logger(__name__)

# ---------- Firestore helpers ----------

_firestore_client = None

def _get_firestore():
    global _firestore_client
    if _firestore_client is not None:
        return _firestore_client
    from google.cloud import firestore  # type: ignore
    project_id = os.getenv("FIRESTORE_PROJECT_ID") or None
    _firestore_client = firestore.Client(project=project_id)
    return _firestore_client

def firestore_patch(job_id: str, patch: Dict[str, Any]) -> None:
    """Best-effort patch. Never kills the job if telemetry fails."""
    try:
        db = _get_firestore()
        ref = db.collection("batch_jobs").document(job_id)
        patch = dict(patch)
        patch["updatedAt"] = int(time.time())
        ref.set(patch, merge=True)
    except Exception as e:
        log.warning("Firestore patch failed job_id=%s err=%s", job_id, e)

# ---------- Domain: per-row analysis ----------

def analyze_one_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """
    TODO: replace with your real single-compound pipeline.
    Current: minimal validation only.
    """
    smiles = row.get("smiles")
    if not smiles or not isinstance(smiles, str):
        raise ValueError("Row missing 'smiles'")
    # Here you should validate with RDKit + compute descriptors
    time.sleep(0.05)
    return {"smiles": smiles, "ok": True}

# ---------- Task ----------

@celery_app.task(name="app.chembind.tasks.process_batch_job", bind=True, max_retries=0)
def process_batch_job(self, payload: Dict[str, Any]) -> Dict[str, Any]:
    job_id = payload.get("job_id")
    rows: List[Dict[str, Any]] = payload.get("rows", [])
    attempt_cap = int(payload.get("attempt_cap", 3))
    attempt = int(payload.get("attempt", 1))

    if not job_id or not isinstance(job_id, str):
        raise ValueError("payload.job_id is required (string)")
    if not isinstance(rows, list):
        raise ValueError("payload.rows must be a list")

    total = len(rows)
    processed = succeeded = failed = 0
    results: List[Dict[str, Any]] = []
    errors: List[Dict[str, Any]] = []

    firestore_patch(job_id, {
        "status": "running",
        "attempt": attempt,
        "attemptCap": attempt_cap,
        "total": total,
        "processed": 0,
        "succeeded": 0,
        "failed": 0,
        "progress": 0.0,
        "startedAt": payload.get("startedAt") or int(time.time()),
    })

    try:
        progress_every = int(os.getenv("BATCH_PROGRESS_EVERY", "10"))

        for idx, row in enumerate(rows):
            processed += 1
            try:
                out = analyze_one_row(row)
                results.append({"index": idx, "result": out})
                succeeded += 1
            except Exception as e:
                failed += 1
                errors.append({"index": idx, "error": str(e)})

            if processed % progress_every == 0 or processed == total:
                firestore_patch(job_id, {
                    "processed": processed,
                    "succeeded": succeeded,
                    "failed": failed,
                    "progress": processed / max(1, total),
                })

        firestore_patch(job_id, {
            "status": "completed",
            "processed": processed,
            "succeeded": succeeded,
            "failed": failed,
            "progress": 1.0,
            "completedAt": int(time.time()),
        })

        return {
            "job_id": job_id,
            "status": "completed",
            "processed": processed,
            "succeeded": succeeded,
            "failed": failed,
            "results": results,
            "errors": errors,
        }

    except SoftTimeLimitExceeded:
        firestore_patch(job_id, {
            "status": "timed_out",
            "processed": processed,
            "succeeded": succeeded,
            "failed": failed,
            "progress": processed / max(1, total),
            "error": "SoftTimeLimitExceeded",
        })
        return {
            "job_id": job_id,
            "status": "timed_out",
            "processed": processed,
            "succeeded": succeeded,
            "failed": failed,
            "results": results,
            "errors": errors,
        }

    except Exception as e:
        firestore_patch(job_id, {
            "status": "failed",
            "processed": processed,
            "succeeded": succeeded,
            "failed": failed,
            "error": str(e),
        })

        # Attempt cap: requeue the same payload (you can optimize to requeue remaining rows later)
        if attempt < attempt_cap:
            new_payload = dict(payload)
            new_payload["attempt"] = attempt + 1
            celery_app.send_task("app.chembind.tasks.process_batch_job", args=[new_payload])
            firestore_patch(job_id, {"status": "requeued", "attempt": attempt + 1})

        raise
