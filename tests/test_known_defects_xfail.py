from __future__ import annotations

from datetime import datetime, timezone

import pytest

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import Visibilities


@pytest.mark.xfail(strict=True, reason="Known defect: negative custom angles do not currently raise ValueError.")
def test_set_custom_prayer_angles_negative_value_should_raise(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(ValueError, match="must be greater than 0"):
        location.set_custom_prayer_angles(fajr_angle=-5.0)


@pytest.mark.xfail(strict=True, reason="Known defect: visibilities(days=...) type error message path is broken.")
def test_visibilities_bad_days_type_should_raise_clear_typeerror(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(TypeError, match="'days' must be of type `int`"):
        location.visibilities(days="3", criterion=1)  # type: ignore[arg-type]


@pytest.mark.xfail(strict=True, reason="Known defect: Visibilities.__str__ returns early when q == 0.")
def test_visibilities_str_with_zero_q_should_include_context_lines() -> None:
    vis = Visibilities(
        criterion="Yallop",
        dates=(datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),),
        q_values=(0.0,),
        classifications=("A: Easily visible.",),
    )
    rendered = str(vis)
    assert "Visibility of New Moon Crescent" in rendered
    assert "Criterion" in rendered
