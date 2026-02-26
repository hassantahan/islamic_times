from __future__ import annotations

import math

from islamic_times import calculation_equations as ce
from islamic_times.it_dataclasses import Angle, Distance, DistanceUnits, RightAscension


def test_sin_and_cos_degree_helpers() -> None:
    assert math.isclose(float(ce.sin(30)), 0.5, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(float(ce.cos(60)), 0.5, rel_tol=0.0, abs_tol=1e-12)


def test_calculate_angle_diff_identity_and_opposite_points() -> None:
    assert math.isclose(ce.calculate_angle_diff(10, 15, 10, 15), 0.0, abs_tol=1e-9)
    assert math.isclose(ce.calculate_angle_diff(0, 0, 180, 0), 180.0, abs_tol=1e-9)


def test_haversine_same_location_gives_zero_distance() -> None:
    distance_km, bearing_deg = ce.haversine(43.65107, -79.347015, 43.65107, -79.347015)
    assert math.isclose(distance_km, 0.0, abs_tol=1e-9)
    assert isinstance(bearing_deg, float)


def test_get_cardinal_direction_boundaries() -> None:
    assert ce.get_cardinal_direction(0) == "N"
    assert ce.get_cardinal_direction(22.5) == "NNE"
    assert ce.get_cardinal_direction(90) == "E"
    assert ce.get_cardinal_direction(225) == "SW"


def test_angle_diff_signed_normalization() -> None:
    assert math.isclose(ce.angle_diff(350, 10), 20.0, abs_tol=1e-12)
    assert math.isclose(ce.angle_diff(10, 350), -20.0, abs_tol=1e-12)
    assert math.isclose(ce.angle_diff(10, 190), 180.0, abs_tol=1e-12)


def test_interpolation_wraps_into_0_360_range() -> None:
    interpolated = ce.interpolation(0.5, 350, 10, 20)
    assert 0.0 <= interpolated < 360.0


def test_correct_ra_dec_returns_ra_and_dec_types() -> None:
    ascension, declination = ce.correct_ra_dec(
        ra=RightAscension(5.5),
        dec=Angle(20),
        lha=Angle(30),
        parallax=Angle(1),
        lat=Angle(43),
        elev=Distance(100, DistanceUnits.METRE),
    )
    assert isinstance(ascension, RightAscension)
    assert isinstance(declination, Angle)


def test_geocentric_horizontal_coordinates_output_ranges() -> None:
    altitude, azimuth = ce.geocentric_horizontal_coordinates(
        observer_latitude=Angle(43.65107),
        body_declination=Angle(15),
        body_lha=Angle(30),
    )
    assert isinstance(altitude, Angle)
    assert isinstance(azimuth, Angle)
    assert -90.0 <= altitude.decimal <= 90.0
    assert 0.0 <= azimuth.decimal < 360.0
