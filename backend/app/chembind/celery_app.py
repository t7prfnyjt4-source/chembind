# app/chembind/celery_app.py
from __future__ import annotations

import os
from celery import Celery


def make_celery() -> Celery:
    # Prefer REDIS_URL, fall back to CELERY_BROKER_URL
    broker = os.getenv("REDIS_URL") or os.getenv("CELERY_BROKER_URL")
    if not broker:
        raise RuntimeError("REDIS_URL is required (e.g. redis://localhost:6379/0)")

    backend = os.getenv("CELERY_RESULT_BACKEND", broker)

    celery_app = Celery(
        "chembind",
        broker=broker,
        backend=backend,
        include=["app.chembind.tasks"],
    )

    # Optional tuning from env
    visibility_timeout = int(os.getenv("CELERY_VISIBILITY_TIMEOUT_SECONDS", "3600"))
    soft_time_limit = int(os.getenv("CELERY_SOFT_TIME_LIMIT_SECONDS", "540"))
    hard_time_limit = int(os.getenv("CELERY_HARD_TIME_LIMIT_SECONDS", "600"))
    prefetch = int(os.getenv("CELERY_PREFETCH_MULTIPLIER", "1"))

    celery_app.conf.update(
        broker_transport_options={"visibility_timeout": visibility_timeout},
        task_soft_time_limit=soft_time_limit,
        task_time_limit=hard_time_limit,
        worker_prefetch_multiplier=prefetch,
        task_acks_late=True,
    )

    return celery_app


celery_app = make_celery()
