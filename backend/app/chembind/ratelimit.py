from __future__ import annotations

import os
import time
from collections import deque
from typing import Callable, Deque, Dict

from fastapi import HTTPException, Request

# Simple in-memory IP limiter (single-instance MVP)
# Env:
#   RATE_LIMIT_RPM=60
#   RATE_LIMIT_BURST=10
RPM = int(os.getenv("RATE_LIMIT_RPM", "60"))
BURST = int(os.getenv("RATE_LIMIT_BURST", "10"))
WINDOW_S = 60.0

_hits: Dict[str, Deque[float]] = {}

def _client_ip(request: Request) -> str:
    return request.client.host if request.client else "unknown"

def build_limiter() -> Callable[[Request], None]:
    limit = RPM + BURST

    async def limiter(request: Request) -> None:
        now = time.time()
        ip = _client_ip(request)

        q = _hits.get(ip)
        if q is None:
            q = deque()
            _hits[ip] = q

        cutoff = now - WINDOW_S
        while q and q[0] < cutoff:
            q.popleft()

        if len(q) >= limit:
            raise HTTPException(status_code=429, detail="Rate limit exceeded")

        q.append(now)

    return limiter
