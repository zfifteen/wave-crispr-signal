from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Protocol, Sequence


@dataclass(frozen=True)
class ComparatorRecord:
    guide_seq: str
    target_context: str


@dataclass(frozen=True)
class ComparatorResult:
    ok: bool
    message: str
    details: Dict[str, object]


class ComparatorAdapter(Protocol):
    name: str
    version: str

    def self_check(self) -> ComparatorResult:
        ...

    def provenance(self) -> Dict[str, object]:
        ...

    def predict_batch(self, records: Sequence[ComparatorRecord]) -> List[float]:
        ...


def resolve_gate_root() -> Path:
    return Path(__file__).resolve().parents[1]
