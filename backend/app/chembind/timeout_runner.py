from __future__ import annotations

import multiprocessing as mp
from dataclasses import dataclass
from typing import Any, Callable, Dict, Tuple


class TimeoutError(RuntimeError):
    pass


class WorkerError(RuntimeError):
    pass


@dataclass(frozen=True)
class TimeoutConfig:
    seconds: float = 2.5


def _worker(fn: Callable[..., Any], args: Tuple[Any, ...], kwargs: Dict[str, Any], q: mp.Queue) -> None:
    try:
        q.put(("ok", fn(*args, **kwargs)))
    except Exception as e:
        q.put(("exc", e.__class__.__name__, str(e)))


def run_with_timeout(
    fn: Callable[..., Any],
    cfg: TimeoutConfig,
    *args: Any,
    **kwargs: Any,
) -> Any:
    """
    IMPORTANT: fn MUST be a top-level function (picklable).
    Args/kwargs must also be picklable.
    """
    ctx = mp.get_context("spawn")
    q: mp.Queue = ctx.Queue()

    p = ctx.Process(target=_worker, args=(fn, args, kwargs, q), daemon=True)
    p.start()
    p.join(timeout=cfg.seconds)

    if p.is_alive():
        p.terminate()
        p.join(timeout=0.2)
        raise TimeoutError(f"Computation exceeded {cfg.seconds}s")

    if q.empty():
        raise WorkerError("Worker exited without returning a result")

    status = q.get()
    if status[0] == "ok":
        return status[1]

    # ("exc", exc_name, exc_msg)
    _, exc_name, exc_msg = status

    # bubble up known validation errors so FastAPI handler returns INVALID_SMILES
    if exc_name == "SmilesValidationError":
        # import here to avoid circular imports
        from app.chembind.rdkit_safe import SmilesValidationError
        raise SmilesValidationError(exc_msg)

    raise WorkerError(f"{exc_name}: {exc_msg}")
