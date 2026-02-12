from __future__ import annotations

from datetime import datetime, timezone

import pytest

import islamic_times.astro_core as fast_astro
from islamic_times.islamic_times import ITLocation


@pytest.mark.parametrize("longitude", [-179.0, -178.5, -178.0])
def test_public_sun_transit_stays_on_observer_date_at_dateline(longitude: float) -> None:
    """Sun transit for the reference date should not roll into the next UTC date."""
    location = ITLocation(
        latitude=0.0,
        longitude=longitude,
        elevation=10.0,
        temperature=15.0,
        pressure=101.325,
        date=datetime(2024, 1, 8, 12, 0, tzinfo=timezone.utc),
        method="ISNA",
        find_local_tz=False,
    )
    observer_date = location.dates_times().date
    transit = location.sun().sun_transit

    assert (transit.date() - observer_date.date()).days == 0


@pytest.mark.parametrize("longitude", [-179.0, -178.5, -178.0])
def test_public_moon_transit_stays_on_observer_date_at_dateline(longitude: float) -> None:
    """Moon transit for the reference date should not roll into the next UTC date."""
    location = ITLocation(
        latitude=0.0,
        longitude=longitude,
        elevation=10.0,
        temperature=15.0,
        pressure=101.325,
        date=datetime(2024, 1, 8, 12, 0, tzinfo=timezone.utc),
        method="ISNA",
        find_local_tz=False,
    )
    observer_date = location.dates_times().date
    transit = location.moon().moon_transit

    assert (transit.date() - observer_date.date()).days == 0


def test_native_sun_transit_not_next_day_for_dateline_case() -> None:
    """Regression for old +1 day behavior near longitude -179."""
    input_dt = datetime(2024, 1, 8, 12, 0, tzinfo=timezone.utc)
    jd = fast_astro.gregorian_to_jd(input_dt, 0.0)
    delta_t = fast_astro.delta_t_approx(input_dt.year, input_dt.month)

    transit = fast_astro.find_sun_transit(
        jd,
        delta_t,
        0.0,
        -179.0,
        10.0,
        15.0,
        101.325,
        0.0,
    )

    assert transit.date() == input_dt.date()


def test_native_moon_transit_not_next_day_for_dateline_case() -> None:
    """Regression for old +1 day behavior near longitude -179."""
    input_dt = datetime(2024, 1, 8, 12, 0, tzinfo=timezone.utc)
    jd = fast_astro.gregorian_to_jd(input_dt, 0.0)
    delta_t = fast_astro.delta_t_approx(input_dt.year, input_dt.month)
    sentinel = [-123456.0, -123456.0, -123456.0]

    transit = fast_astro.find_moon_transit(
        jd,
        delta_t,
        0.0,
        -179.0,
        10.0,
        15.0,
        101.325,
        0.0,
        sentinel,
        sentinel,
    )

    assert transit.date() == input_dt.date()
