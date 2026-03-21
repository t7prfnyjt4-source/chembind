# app/chembind/tier_gating.py
"""User tier gating based on Firebase Custom Claims."""
from __future__ import annotations

from typing import Any, Dict

from fastapi import HTTPException

TIER_HIERARCHY: Dict[str, int] = {
    "basic": 0,
    "pro": 1,
    "professional": 2,
    "team": 3,
    "enterprise": 4,
}


def get_user_tier(user: Dict[str, Any]) -> str:
    """Extract tier from user claims, default to 'basic'."""
    claims = user.get("claims", {})
    return claims.get("tier", "basic")


def check_tier(user: Dict[str, Any], required_tier: str) -> None:
    """
    Raise 402 if user's tier is below required_tier.
    """
    user_tier = get_user_tier(user)
    user_level = TIER_HIERARCHY.get(user_tier, 0)
    required_level = TIER_HIERARCHY.get(required_tier, 0)

    if user_level < required_level:
        raise HTTPException(
            status_code=402,
            detail=f"Requires {required_tier} tier or above (current: {user_tier})",
        )
