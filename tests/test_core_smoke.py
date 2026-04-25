from __future__ import annotations

from datetime import datetime, timedelta, timezone
import math

import pytest

from islamic_times import sun_equations as se
from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import Angle, ObserverInfo, PrayerTimes, Visibilities


def test_public_api_smoke() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=10.0,
        temperature=15.0,
        pressure=101.325,
        date=datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),
        method="ISNA",
        find_local_tz=False,
    )

    observer = location.observer()
    prayers = location.prayer_times()

    assert isinstance(observer, ObserverInfo)
    assert isinstance(prayers, PrayerTimes)
    assert prayers.fajr.name.lower() == "fajr"
    assert prayers.sunrise.time is not None


def test_visibilities_smoke() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=10.0,
        date=datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),
        find_local_tz=False,
    )
    visibilities = location.visibilities(days=1, criterion=1)

    assert isinstance(visibilities, Visibilities)
    assert len(visibilities.dates) == 1
    assert len(visibilities.q_values) == 1
    assert len(visibilities.classifications) == 1


def test_visibilities_smoke_shaukat() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=10.0,
        date=datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),
        find_local_tz=False,
    )
    visibilities = location.visibilities(days=1, criterion=2)

    assert isinstance(visibilities, Visibilities)
    assert visibilities.criterion == "Shaukat"
    assert len(visibilities.dates) == 1
    assert len(visibilities.q_values) == 1
    assert len(visibilities.classifications) == 1


def test_itlocation_uses_modern_timescales_for_toronto_reference_case() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=150.0,
        temperature=10.0,
        pressure=101.325,
        date=datetime(2026, 3, 19, 19, 58, 0, tzinfo=timezone(timedelta(hours=-5))),
        method="ISNA",
        find_local_tz=False,
    )

    date_info = location.dates_times()

    assert date_info.deltaT == pytest.approx(69.0928431609, abs=0.02)
    assert date_info.tt_utc_seconds == pytest.approx(69.184, abs=1e-6)
    assert date_info.ut1_utc_seconds == pytest.approx(0.0911568391, abs=0.02)


def _bennett_refraction_lift(altitude_deg: float) -> float:
    return 0.017 / math.tan(math.radians(altitude_deg + 10.3 / (altitude_deg + 5.11)))


def test_sun_refraction_clamps_below_minus_one_degree_and_suninfo_uses_apparent() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=150.0,
        temperature=10.0,
        pressure=101.325,
        date=datetime(2026, 4, 6, 4, 0, 0, tzinfo=timezone.utc),
        method="ISNA",
        find_local_tz=False,
    )

    true_altitude = location.sun_params.true_altitude.decimal
    apparent_altitude = location.sun_params.apparent_altitude.decimal

    assert true_altitude < -1.0
    assert apparent_altitude - true_altitude == pytest.approx(_bennett_refraction_lift(-1.0), abs=1e-9)
    assert location.sun().apparent_altitude.decimal == pytest.approx(apparent_altitude, abs=1e-12)


def test_sun_refraction_uses_true_altitude_above_minus_one_degree() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=150.0,
        temperature=10.0,
        pressure=101.325,
        date=datetime(2026, 4, 6, 10, 50, 0, tzinfo=timezone.utc),
        method="ISNA",
        find_local_tz=False,
    )

    true_altitude = location.sun_params.true_altitude.decimal
    apparent_altitude = location.sun_params.apparent_altitude.decimal

    assert true_altitude > -1.0
    assert apparent_altitude - true_altitude == pytest.approx(_bennett_refraction_lift(true_altitude), abs=1e-9)


def test_custom_solar_angle_event_uses_geometric_altitude_target() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=150.0,
        temperature=10.0,
        pressure=101.325,
        date=datetime(2026, 4, 6, 12, 0, 0, tzinfo=timezone.utc),
        method="ISNA",
        find_local_tz=False,
    )

    event_dt = se.find_proper_suntime(
        location.observer_dateinfo,
        location.observer_info,
        "rise",
        Angle(12.0),
    )
    event_location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=150.0,
        temperature=10.0,
        pressure=101.325,
        date=event_dt,
        method="ISNA",
        find_local_tz=False,
    )

    assert event_location.sun_params.true_altitude.decimal == pytest.approx(-12.0, abs=0.03)
    assert event_location.sun_params.apparent_altitude.decimal > event_location.sun_params.true_altitude.decimal


@pytest.mark.parametrize(
    ("latitude", "longitude", "elevation", "observer_dt"),
    [
        (43.651070, -79.347015, 150.0, datetime(2026, 4, 6, 4, 0, 0, tzinfo=timezone.utc)),
        (51.507222, -0.127500, 11.0, datetime(2026, 4, 6, 6, 58, 0, tzinfo=timezone.utc)),
        (21.389082, 39.857910, 277.0, datetime(2026, 4, 6, 19, 30, 0, tzinfo=timezone.utc)),
    ],
)
def test_moon_refraction_clamps_below_minus_one_degree(
    latitude: float,
    longitude: float,
    elevation: float,
    observer_dt: datetime,
) -> None:
    location = ITLocation(
        latitude=latitude,
        longitude=longitude,
        elevation=elevation,
        temperature=10.0,
        pressure=101.325,
        date=observer_dt,
        method="ISNA",
        find_local_tz=False,
    )

    true_altitude = location.moon_params.true_altitude.decimal
    apparent_altitude = location.moon().apparent_altitude.decimal

    assert true_altitude < -1.0
    assert apparent_altitude - true_altitude == pytest.approx(_bennett_refraction_lift(-1.0), abs=1e-9)


def test_moon_refraction_uses_true_altitude_above_minus_one_degree() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=150.0,
        temperature=10.0,
        pressure=101.325,
        date=datetime(2026, 3, 19, 23, 58, 0, tzinfo=timezone.utc),
        method="ISNA",
        find_local_tz=False,
    )

    true_altitude = location.moon_params.true_altitude.decimal
    apparent_altitude = location.moon().apparent_altitude.decimal

    assert true_altitude > -1.0
    assert apparent_altitude - true_altitude == pytest.approx(_bennett_refraction_lift(true_altitude), abs=1e-9)
