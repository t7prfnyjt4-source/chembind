#!/usr/bin/env bash
set -euo pipefail

# Only load /app/.env for keys NOT already set by docker-compose env.
dotenv="/app/.env"
if [ -f "$dotenv" ]; then
  while IFS='=' read -r key val; do
    [ -z "${key:-}" ] && continue
    case "$key" in \#*) continue ;; esac

    val="${val%$'\r'}"
    val="${val%\"}"; val="${val#\"}"
    val="${val%\'}"; val="${val#\'}"

    if [ -z "${!key-}" ]; then
      export "$key=$val"
    fi
  done < <(grep -v '^[[:space:]]*$' "$dotenv" | grep -v '^[[:space:]]*#')
fi

: "${CELERY_APP:=app.chembind.celery_app:celery_app}"
: "${CELERY_CONCURRENCY:=2}"
: "${CELERY_MAX_TASKS_PER_CHILD:=50}"
: "${CELERY_LOGLEVEL:=INFO}"
: "${CELERY_PREFETCH_MULTIPLIER:=1}"
: "${CELERY_POOL:=prefork}"
: "${CELERY_QUEUES:=default}"

exec micromamba run -n app celery -A "$CELERY_APP" worker \
  --loglevel="$CELERY_LOGLEVEL" \
  --concurrency="$CELERY_CONCURRENCY" \
  --max-tasks-per-child="$CELERY_MAX_TASKS_PER_CHILD" \
  --prefetch-multiplier="$CELERY_PREFETCH_MULTIPLIER" \
  --pool="$CELERY_POOL" \
  -Q "$CELERY_QUEUES"

