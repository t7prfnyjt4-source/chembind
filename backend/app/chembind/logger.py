from __future__ import annotations

import contextvars
import logging
import sys
import uuid
from typing import Optional

from pythonjsonlogger import jsonlogger

# Async-safe request id storage (per request)
request_id_ctx_var: contextvars.ContextVar[Optional[str]] = contextvars.ContextVar(
    "request_id", default=None
)

class RequestIdFilter(logging.Filter):
    """
    Inject a request_id attribute into every log record.
    If no request id exists, it becomes "-".
    """
    def filter(self, record: logging.LogRecord) -> bool:
        record.request_id = request_id_ctx_var.get() or "-"
        return True

def set_request_id(request_id: Optional[str]) -> None:
    request_id_ctx_var.set(request_id)

def new_request_id() -> str:
    return str(uuid.uuid4())

def setup_json_logging(level: str = "INFO") -> None:
    """
    Configure root logger to emit JSON to stdout.
    Adds request_id to every record via RequestIdFilter.
    """
    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(level.upper())

    handler = logging.StreamHandler(sys.stdout)
    formatter = jsonlogger.JsonFormatter(
        "%(asctime)s %(levelname)s %(name)s %(message)s %(request_id)s"
    )
    handler.setFormatter(formatter)
    handler.addFilter(RequestIdFilter())

    root.addHandler(handler)

# Backwards-compat alias (Segment 4 middleware expects request_id_var)
request_id_var = request_id_ctx_var

# Backwards-compat alias
request_id_var = request_id_ctx_var

# Backwards-compat alias
request_id_var = request_id_ctx_var
