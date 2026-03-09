from __future__ import annotations

from typing import Dict

from .base import ComparatorAdapter
from .crisprpred.adapter import CRISPRpredAdapter


_RESERVED_SLOTS = {
    "baseline_c": "crisprpred",
    "baseline_d": None,
}


def required_comparator_slots() -> Dict[str, str | None]:
    return dict(_RESERVED_SLOTS)


def build_required_comparators() -> Dict[str, ComparatorAdapter]:
    return {
        "baseline_c": CRISPRpredAdapter(),
    }
