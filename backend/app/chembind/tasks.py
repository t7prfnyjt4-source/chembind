# app/chembind/tasks.py
from __future__ import annotations

from typing import Any, Dict
import os

from app.chembind.rdkit_safe import compute_descriptors, SmilesValidationError, RdkitLimits
from app.chembind.firestore_repo import FirestoreRepo
from app.chembind.celery_app import celery_app


@celery_app.task(name="chembind.process_batch_job")
def process_batch_job(uid: str, job_id: str) -> Dict[str, Any]:
    """
    Process a batch job:
    - load input rows from Firestore (users/{uid}/jobs/{jobId}/items)
    - update job status running
    - compute descriptors for each row
    - write per-row item results
    - update job totals + finished status
    """
    repo = FirestoreRepo()

    # ✅ Load input rows for this job from Firestore
    rows = repo.list_job_items(uid, job_id)

    max_smiles_len = int(os.getenv("MAX_SMILES_LEN", "500"))
    max_atoms = int(os.getenv("MAX_ATOMS", "200"))
    limits = RdkitLimits(max_smiles_len=max_smiles_len, max_atoms=max_atoms)

    total = len(rows)
    processed = 0
    success = 0
    failed = 0

    # mark running
    repo.update_job(
        uid,
        job_id,
        {
            "status": "running",
            "total": total,
            "processed": 0,
            "successCount": 0,
            "failureCount": 0,
        },
    )

    for idx, row in enumerate(rows):
        processed += 1
        smiles = (row.get("smiles") or "").strip()

        try:
            if not smiles:
                raise ValueError("Missing smiles")

            desc = compute_descriptors(smiles, limits=limits)
            success += 1

            repo.write_job_item(
                uid,
                job_id,
                item_id=str(idx),
                payload={
                    "index": idx,
                    "smiles": smiles,
                    "ok": True,
                    "descriptors": desc,
                },
            )

        except (SmilesValidationError, ValueError) as e:
            failed += 1
            repo.write_job_item(
                uid,
                job_id,
                item_id=str(idx),
                payload={
                    "index": idx,
                    "smiles": smiles,
                    "ok": False,
                    "error": str(e),
                },
            )

        except Exception as e:
            failed += 1
            repo.write_job_item(
                uid,
                job_id,
                item_id=str(idx),
                payload={
                    "index": idx,
                    "smiles": smiles,
                    "ok": False,
                    "error": f"INTERNAL: {repr(e)}",
                },
            )

        # progress updates
        every = int(os.getenv("BATCH_PROGRESS_EVERY", "10"))
        if processed % every == 0 or processed == total:
            repo.update_job(
                uid,
                job_id,
                {
                    "processed": processed,
                    "successCount": success,
                    "failureCount": failed,
                },
            )

    repo.update_job(
        uid,
        job_id,
        {
            "status": "finished",
            "processed": processed,
            "successCount": success,
            "failureCount": failed,
        },
    )

    return {
        "jobId": job_id,
        "status": "finished",
        "total": total,
        "processed": processed,
        "success": success,
        "failed": failed,
    }
