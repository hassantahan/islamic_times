from __future__ import annotations

import pytest

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import MeccaInfo, MoonInfo, SunInfo, Visibilities


def test_set_prayer_method_is_case_insensitive(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    location.set_prayer_method("mwl", asr_type=0)
    assert "Muslim World League" in location.prayer_times().method.name


def test_set_prayer_method_rejects_invalid_asr_type(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(ValueError, match="'asr_type' must be either 0 or 1"):
        location.set_prayer_method("ISNA", asr_type=2)


def test_set_custom_prayer_angles_updates_method_to_custom(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    location.set_custom_prayer_angles(fajr_angle=17.0, maghrib_angle=4.0, isha_angle=15.0)

    method = location.prayer_times().method
    assert method.name == "Custom"
    assert method.fajr_angle.decimal == pytest.approx(17.0)
    assert method.maghrib_angle.decimal == pytest.approx(4.0)
    assert method.isha_angle.decimal == pytest.approx(15.0)


def test_set_custom_prayer_angles_rejects_non_numeric_input(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(ValueError, match="must be a number"):
        location.set_custom_prayer_angles(fajr_angle="bad")  # type: ignore[arg-type]


def test_set_custom_prayer_angles_negative_value_should_raise(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(ValueError, match="must be greater than 0"):
        location.set_custom_prayer_angles(fajr_angle=-5.0)


def test_set_asr_type_updates_method(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    location.set_asr_type(1)
    assert location.prayer_times().method.asr_type == 1
    assert location.prayer_times().method.name == "Custom"


def test_set_asr_type_rejects_invalid_value(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(ValueError, match="'asr_type' must be either 0 or 1"):
        location.set_asr_type(9)


def test_set_midnight_type_updates_method(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    location.set_midnight_type(1)
    assert location.prayer_times().method.midnight_type == 1
    assert location.prayer_times().method.name == "Custom"


def test_set_midnight_type_rejects_invalid_value(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(ValueError, match="'midnight_type' must be either 0 or 1"):
        location.set_midnight_type(9)


def test_set_extreme_latitude_rule_validates_input(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)

    with pytest.raises(TypeError, match="'rule' must be a string type"):
        location.set_extreme_latitude_rule(7)  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="is not a valid rule"):
        location.set_extreme_latitude_rule("INVALID")

    location.set_extreme_latitude_rule("ANGLEBASED")
    assert location.prayer_times().method.extreme_lats == "ANGLEBASED"


def test_mecca_returns_typed_result(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    mecca = location.mecca()

    assert isinstance(mecca, MeccaInfo)
    assert mecca.distance.value > 0
    assert 0.0 <= mecca.angle.decimal < 360.0
    assert isinstance(mecca.cardinal, str)


def test_sun_and_moon_accessors_return_typed_results(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    assert isinstance(location.sun(), SunInfo)
    assert isinstance(location.moon(), MoonInfo)


def test_moonphases_returns_four_named_events(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    phases = location.moonphases()

    assert len(phases) == 4
    assert [name for name, _ in phases] == ["New Moon", "First Quarter", "Full Moon", "Last Quarter"]


def test_visibilities_valid_small_input(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    visibilities = location.visibilities(days=1, criterion=1)

    assert isinstance(visibilities, Visibilities)
    assert len(visibilities.dates) == 1


def test_visibilities_rejects_invalid_days_and_criterion(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)

    with pytest.raises(ValueError, match="'days' must be greater than 0"):
        location.visibilities(days=0, criterion=1)

    with pytest.raises(ValueError, match="'criterion' must be either 0 or 1"):
        location.visibilities(days=1, criterion=3)

    with pytest.raises(TypeError, match="'criterion' must be of type `int`"):
        location.visibilities(days=1, criterion="1")  # type: ignore[arg-type]


def test_visibilities_rejects_non_integer_days(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(TypeError, match="'days' must be of type `int`"):
        location.visibilities(days="3", criterion=1)  # type: ignore[arg-type]
