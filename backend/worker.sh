#!/usr/bin/env bash
set -euo pipefail

: "${CELERY_CONCURRENCY:=2}"
: "${CELERY_MAX_TASKS_PER_CHILD:=50}"
: "${CELERY_LOGLEVEL:=INFO}"
: "${CELERY_PREFETCH_MULTIPLIER:=1}"

export CELERYD_PREFETCH_MULTIPLIER="$CELERY_PREFETCH_MULTIPLIER"

exec celery \
  -A app.chembind.celery_worker \
  worker \
  --loglevel="$CELERY_LOGLEVEL" \
  --concurrency="$CELERY_CONCURRENCY" \
  --max-tasks-per-child="$CELERY_MAX_TASKS_PER_CHILD"
