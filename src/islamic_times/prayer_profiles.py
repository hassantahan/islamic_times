"""Reusable mosque prayer methods and their local timetable offsets."""

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from datetime import timedelta
from math import isfinite
from types import MappingProxyType
from typing import Any

from islamic_times.it_dataclasses import PRAYER_KEYS, PUBLIC_SCHEMA_VERSION, Angle, Prayer, PrayerMethod, PrayerTimes
from islamic_times.prayer_times import HIGH_LATITUDE_RULES, get_prayer_method
from islamic_times.time_equations import add_elapsed_time


@dataclass(frozen=True, slots=True)
class MosqueProfile:
    """A named calculation method with local adjustments to prayer times.

    Parameters
    ----------
    name : str
        Display name for the mosque or timetable.
    method : PrayerMethod
        Calculation settings to reuse.
    offset_minutes : Mapping[str, int], optional
        Signed whole-minute adjustments keyed by ``fajr``, ``zuhr``,
        ``asr``, ``maghrib``, or ``isha``. Omitted keys have no adjustment.
    """

    name: str
    method: PrayerMethod
    offset_minutes: Mapping[str, int] = field(default_factory=dict, hash=False)

    def __post_init__(self) -> None:
        """Validate and snapshot the caller's offsets."""
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("'name' must be a non-empty string.")
        if not isinstance(self.method, PrayerMethod):
            raise TypeError("'method' must be a PrayerMethod.")
        if not isinstance(self.offset_minutes, Mapping):
            raise TypeError("'offset_minutes' must be a mapping.")

        offsets = dict(self.offset_minutes)
        for key, offset in offsets.items():
            if key not in PRAYER_KEYS:
                raise ValueError(f"Unknown prayer key {key!r}. Use one of: {', '.join(PRAYER_KEYS)}.")
            if type(offset) is not int:
                raise TypeError(f"The offset for {key!r} must be a whole number of minutes.")
            try:
                timedelta(minutes=offset)
            except OverflowError:
                raise ValueError(f"The offset for {key!r} is too large.") from None
        object.__setattr__(self, "offset_minutes", MappingProxyType(offsets))

    @classmethod
    def from_preset(
        cls, name: str, method_key: str, *, offset_minutes: Mapping[str, int] | None = None,
    ) -> "MosqueProfile":
        """Build a profile from a built-in calculation method or its alias."""
        offsets = {} if offset_minutes is None else offset_minutes
        return cls(name=name, method=get_prayer_method(method_key), offset_minutes=offsets)

    def to_dict(self) -> dict[str, Any]:
        """Export the complete profile to a JSON-compatible dictionary."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "MosqueProfile",
            "name": self.name,
            "method": self.method.to_dict(),
            "offset_minutes": dict(self.offset_minutes),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "MosqueProfile":
        """Restore a profile produced by :meth:`to_dict`, validating its settings."""
        if not isinstance(payload, Mapping):
            raise TypeError("The mosque profile must be a mapping.")
        if payload.get("type") != "MosqueProfile" or payload.get("schema_version") != PUBLIC_SCHEMA_VERSION:
            raise ValueError("Unsupported mosque profile format or schema version.")
        return cls(
            name=payload.get("name"),
            method=_method_from_dict(payload.get("method")),
            offset_minutes=payload.get("offset_minutes", {}),
        )

    def _apply_offsets(self, prayers: PrayerTimes) -> PrayerTimes:
        """Apply offsets to available prayer times after calculation."""
        adjusted: dict[str, Prayer] = {}
        for key, offset in self.offset_minutes.items():
            prayer: Prayer = getattr(prayers, key)
            when = prayer.when
            if offset and when is not None:
                adjusted[key] = replace(
                    prayer, time=add_elapsed_time(when, timedelta(minutes=offset)), offset_minutes=offset,
                )
        return replace(prayers, mosque_profile_name=self.name, **adjusted)


def _method_from_dict(payload: object) -> PrayerMethod:
    """Validate the serialized PrayerMethod embedded in a profile."""
    if not isinstance(payload, Mapping):
        raise TypeError("The profile's prayer method must be a mapping.")
    if payload.get("type") != "PrayerMethod" or payload.get("schema_version") != PUBLIC_SCHEMA_VERSION:
        raise ValueError("Unsupported prayer method format or schema version.")

    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("The prayer method name must be a non-empty string.")
    keys = payload.get("keys")
    if keys is not None and (
        not isinstance(keys, (tuple, list)) or any(not isinstance(key, str) for key in keys)
    ):
        raise ValueError("The prayer method keys must be a list of strings or null.")

    angles: dict[str, Angle] = {}
    for key in ("fajr_angle", "maghrib_angle", "isha_angle"):
        degrees = payload.get(f"{key}_degrees")
        if type(degrees) not in (float, int):
            raise ValueError(f"The prayer method {key} must be a finite non-negative number of degrees.")
        try:
            value = float(degrees)
        except OverflowError:
            raise ValueError(f"The prayer method {key} is too large.") from None
        if not isfinite(value) or value < 0:
            raise ValueError(f"The prayer method {key} must be a finite non-negative number of degrees.")
        angles[key] = Angle(value)

    asr_type, midnight_type = payload.get("asr_type"), payload.get("midnight_type")
    if not isinstance(asr_type, int) or asr_type not in (0, 1, 2):
        raise ValueError("Prayer method asr_type must be 0, 1, or 2.")
    if not isinstance(midnight_type, int) or midnight_type not in (0, 1):
        raise ValueError("Prayer method midnight_type must be 0 or 1.")
    extreme_lats = payload.get("extreme_lats")
    if extreme_lats not in HIGH_LATITUDE_RULES:
        raise ValueError("Unknown high-latitude rule in the prayer method.")

    isha_minutes = payload.get("isha_offset_minutes")
    ramadan_minutes = payload.get("ramadan_isha_offset_minutes")
    for offset in (isha_minutes, ramadan_minutes):
        if offset is not None and (type(offset) is not int or offset < 0):
            raise ValueError("Fixed ʿIshāʾ intervals must be non-negative whole minutes or null.")
        if offset is not None:
            try:
                timedelta(minutes=offset)
            except OverflowError:
                raise ValueError("The fixed ʿIshāʾ interval is too large.") from None

    return PrayerMethod(
        name=name,
        keys=tuple(keys) if keys is not None else None,
        asr_type=asr_type,
        midnight_type=midnight_type,
        extreme_lats=extreme_lats,
        isha_offset_minutes=isha_minutes,
        ramadan_isha_offset_minutes=ramadan_minutes,
        **angles,
    )
