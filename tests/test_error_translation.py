from __future__ import annotations

import pytest

import islamic_times.astro_core as fast_astro
from islamic_times import moon_equations as me
from islamic_times import prayer_times as pt
from islamic_times import sun_equations as se
from islamic_times.islamic_times import ITLocation


def test_find_proper_suntime_maps_native_valueerror_to_arithmeticerror(
    monkeypatch: pytest.MonkeyPatch, toronto_observer_kwargs: dict[str, object]
) -> None:
    location = ITLocation(**toronto_observer_kwargs)

    def raise_value_error(*args: object, **kwargs: object) -> None:
        raise ValueError("year -9999 is out of range")

    monkeypatch.setattr(fast_astro, "find_proper_suntime", raise_value_error)
    with pytest.raises(ArithmeticError, match="Sun event does not exist"):
        se.find_proper_suntime(location.observer_dateinfo, location.observer_info, "rise")


def test_find_proper_suntime_does_not_mask_unexpected_runtimeerror(
    monkeypatch: pytest.MonkeyPatch, toronto_observer_kwargs: dict[str, object]
) -> None:
    location = ITLocation(**toronto_observer_kwargs)

    def raise_runtime_error(*args: object, **kwargs: object) -> None:
        raise RuntimeError("unexpected internal failure")

    monkeypatch.setattr(fast_astro, "find_proper_suntime", raise_runtime_error)
    with pytest.raises(RuntimeError, match="unexpected internal failure"):
        se.find_proper_suntime(location.observer_dateinfo, location.observer_info, "set")


def test_find_proper_moontime_maps_native_valueerror_to_arithmeticerror(
    monkeypatch: pytest.MonkeyPatch, toronto_observer_kwargs: dict[str, object]
) -> None:
    location = ITLocation(**toronto_observer_kwargs)

    def raise_value_error(*args: object, **kwargs: object) -> None:
        raise ValueError("year -9999 is out of range")

    monkeypatch.setattr(fast_astro, "find_proper_moontime", raise_value_error)
    with pytest.raises(ArithmeticError, match="Moon event does not exist"):
        me.find_proper_moontime(location.observer_dateinfo, location.observer_info, "set")


def test_find_proper_moontime_does_not_mask_unexpected_runtimeerror(
    monkeypatch: pytest.MonkeyPatch, toronto_observer_kwargs: dict[str, object]
) -> None:
    location = ITLocation(**toronto_observer_kwargs)

    def raise_runtime_error(*args: object, **kwargs: object) -> None:
        raise RuntimeError("unexpected internal failure")

    monkeypatch.setattr(fast_astro, "find_proper_moontime", raise_runtime_error)
    with pytest.raises(RuntimeError, match="unexpected internal failure"):
        me.find_proper_moontime(location.observer_dateinfo, location.observer_info, "set")


def test_calculate_prayer_times_does_not_mask_unexpected_midpoint_error(
    monkeypatch: pytest.MonkeyPatch, toronto_observer_kwargs: dict[str, object]
) -> None:
    location = ITLocation(**toronto_observer_kwargs)

    def raise_runtime_error(*args: object, **kwargs: object) -> None:
        raise RuntimeError("midpoint failed unexpectedly")

    monkeypatch.setattr(pt.te, "time_midpoint", raise_runtime_error)
    with pytest.raises(RuntimeError, match="midpoint failed unexpectedly"):
        pt.calculate_prayer_times(
            location.observer_dateinfo,
            location.observer_info,
            location.sun_info,
            location.method,
        )
