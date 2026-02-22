import os
from celery import Celery

def make_celery() -> Celery:
    redis_url = os.getenv("REDIS_URL")
    if not redis_url:
        raise RuntimeError("REDIS_URL is required (e.g. redis://localhost:6379/0)")

    result_backend = os.getenv("CELERY_RESULT_BACKEND", redis_url)

    celery = Celery(
        "chembind",
        broker=redis_url,
        backend=result_backend,
        include=["app.chembind.tasks"],
    )

    visibility_timeout = int(os.getenv("CELERY_VISIBILITY_TIMEOUT_SECONDS", "3600"))
    soft_limit = int(os.getenv("CELERY_SOFT_TIME_LIMIT_SECONDS", "540"))
    hard_limit = int(os.getenv("CELERY_HARD_TIME_LIMIT_SECONDS", "600"))

    celery.conf.update(
        task_acks_late=True,
        task_reject_on_worker_lost=True,
        worker_prefetch_multiplier=1,
        task_track_started=True,
        task_soft_time_limit=soft_limit,
        task_time_limit=hard_limit,
        broker_transport_options={"visibility_timeout": visibility_timeout},
        task_serializer="json",
        result_serializer="json",
        accept_content=["json"],
        result_expires=86400,
    )

    return celery

celery_app = make_celery()
