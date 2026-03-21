# app/chembind/docking/tasks.py
"""Docking Celery tasks."""
from __future__ import annotations

import logging
import os
import tempfile

from celery.exceptions import SoftTimeLimitExceeded

from app.chembind.celery_app import celery_app
from app.chembind.firestore_repo import FirestoreRepo

logger = logging.getLogger("chembind.docking")

DOCKING_SOFT_TIME_LIMIT = int(os.getenv("DOCKING_SOFT_TIME_LIMIT_S", "120"))


@celery_app.task(
    name="chembind.run_docking_job",
    bind=True,
    soft_time_limit=DOCKING_SOFT_TIME_LIMIT,
    time_limit=DOCKING_SOFT_TIME_LIMIT + 30,
    max_retries=0,
)
def run_docking_job(
    self,
    uid: str,
    job_id: str,
    smiles: str,
    pdb_content: str,
    center: list[float],
    size: list[float],
    exhaustiveness: int = 8,
    num_modes: int = 9,
):
    """Run a docking job: prepare → dock → save results."""
    repo = FirestoreRepo()

    try:
        # Status: preparing
        repo.update_docking_job(uid, job_id, {"status": "preparing"})

        from app.chembind.docking.preparation import (
            prepare_ligand_from_smiles,
            prepare_receptor_from_pdb,
        )
        from app.chembind.docking.engine import run_vina_docking

        with tempfile.TemporaryDirectory() as tmpdir:
            # Prepare ligand
            ligand_pdbqt = prepare_ligand_from_smiles(smiles, tmpdir)

            # Prepare receptor
            receptor_pdbqt = prepare_receptor_from_pdb(pdb_content, tmpdir)

            # Status: docking
            repo.update_docking_job(uid, job_id, {"status": "docking"})

            # Run Vina
            poses = run_vina_docking(
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=ligand_pdbqt,
                center=tuple(center),
                size=tuple(size),
                num_modes=num_modes,
                exhaustiveness=exhaustiveness,
                output_dir=tmpdir,
            )

        # Analyze interactions for top 5 poses (circuit breaker)
        try:
            from app.chembind.docking.interactions import analyze_interactions
            for pose in poses[:5]:
                result = analyze_interactions(pdb_content, pose["atoms"], smiles)
                if result["interactions"]:
                    pose["interactions"] = result["interactions"]
        except Exception as e:
            logger.warning(f"Interaction analysis failed (non-fatal): {e}")

        # Save poses
        repo.save_docking_poses(uid, job_id, poses)

        # Status: completed
        repo.update_docking_job(uid, job_id, {
            "status": "completed",
            "poseCount": len(poses),
            "bestScore": poses[0]["score"] if poses else None,
        })

        logger.info(f"Docking job {job_id} completed: {len(poses)} poses")

    except SoftTimeLimitExceeded:
        logger.warning(f"Docking job {job_id} timed out")
        repo.update_docking_job(uid, job_id, {
            "status": "failed",
            "error": "Docking timed out",
        })

    except Exception as e:
        logger.error(f"Docking job {job_id} failed: {e}")
        repo.update_docking_job(uid, job_id, {
            "status": "failed",
            "error": str(e)[:500],
        })
