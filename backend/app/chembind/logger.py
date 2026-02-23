from __future__ import annotations

import logging
from contextvars import ContextVar

# Async-safe request id storage (per request)
request_id_var: ContextVar[str | None] = ContextVar("request_id", default=None)


class RequestIdFilter(logging.Filter):
    """
    Inject a request_id attribute into every log record.

    This allows log formatters to include %(request_id)s safely.
    If no request id exists, it becomes "-".
    """
    def filter(self, record: logging.LogRecord) -> bool:
        record.request_id = request_id_var.get() or "-"
        return True
