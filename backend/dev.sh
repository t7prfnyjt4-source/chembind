#!/usr/bin/env bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chembind

cd "$(dirname "$0")"
set -a; source .env; set +a

echo "[1/2] Starting worker..."
PYTHONPATH="$(pwd)" ./worker.sh &
WORKER_PID=$!

echo "[2/2] Starting API..."
PYTHONPATH="$(pwd)" uvicorn app.main:app --reload --host 0.0.0.0 --port 8000

trap "kill $WORKER_PID" EXIT
