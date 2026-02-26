from __future__ import annotations

import gc
import tracemalloc

import islamic_times.astro_core as fast_astro


def test_compute_moon_memory_growth_is_bounded() -> None:
    args = (
        2460828.0,
        73.0,
        43.651070,
        -79.347015,
        10.0,
        15.0,
        101.325,
        0.0,
        23.4,
    )

    tracemalloc.start()
    try:
        for _ in range(250):
            fast_astro.compute_moon(*args)

        gc.collect()
        baseline_current, _ = tracemalloc.get_traced_memory()

        for _ in range(5000):
            fast_astro.compute_moon(*args)

        gc.collect()
        current, _ = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert (current - baseline_current) < 250_000
