# backend/app/chembind/celery_app.py
from __future__ import annotations

import os
from celery import Celery


def make_celery() -> Celery:
    """
    Docker-safe Celery configuration.

    Priority order:
    1. CELERY_BROKER_URL / CELERY_RESULT_BACKEND
    2. REDIS_URL
    3. Docker-safe defaults (redis service)
    """

    # --- Resolve broker ---
    broker = (
        os.getenv("CELERY_BROKER_URL")
        or os.getenv("REDIS_URL")
        or "redis://redis:6379/0"
    )

    # --- Resolve backend ---
    backend = os.getenv("CELERY_RESULT_BACKEND")

    if not backend:
        if os.getenv("REDIS_URL"):
            # Use same host but DB 1 for results
            base = os.getenv("REDIS_URL").rsplit("/", 1)[0]
            backend = f"{base}/1"
        else:
            backend = "redis://redis:6379/1"

    celery_app = Celery(
        "chembind",
        broker=broker,
        backend=backend,
        include=["app.chembind.tasks"],
    )

    # --- Optional tuning from environment ---
    celery_app.conf.update(
        broker_transport_options={
            "visibility_timeout": int(
                os.getenv("CELERY_VISIBILITY_TIMEOUT_SECONDS", "3600")
            )
        },
        task_soft_time_limit=int(
            os.getenv("CELERY_SOFT_TIME_LIMIT_SECONDS", "540")
        ),
        task_time_limit=int(
            os.getenv("CELERY_HARD_TIME_LIMIT_SECONDS", "600")
        ),
        worker_prefetch_multiplier=int(
            os.getenv("CELERY_PREFETCH_MULTIPLIER", "1")
        ),
        task_acks_late=True,
        broker_connection_retry_on_startup=True,
    )

    return celery_app


celery_app = make_celery()
