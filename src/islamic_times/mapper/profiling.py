"""Small profiling primitives for mapper pipeline instrumentation."""

from __future__ import annotations

from dataclasses import dataclass, field
from time import perf_counter
from typing import Dict


@dataclass(slots=True)
class StageTimer:
    """Collect wall-clock timings by stage name."""

    _start: Dict[str, float] = field(default_factory=dict)
    elapsed: Dict[str, float] = field(default_factory=dict)

    def start(self, stage: str) -> None:
        self._start[stage] = perf_counter()

    def stop(self, stage: str) -> float:
        t0 = self._start.pop(stage, None)
        if t0 is None:
            return 0.0
        dt = perf_counter() - t0
        self.elapsed[stage] = self.elapsed.get(stage, 0.0) + dt
        return dt

