from __future__ import annotations

from datetime import datetime

from islamic_times.islamic_times import ITLocation


def main() -> None:
    """Small end-to-end example for local verification."""
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=10.0,
        temperature=15.0,
        pressure=101.325,
        date=datetime(2025, 6, 1, 12, 0, 0),
        method="ISNA",
        find_local_tz=True,
    )
    print(location.prayer_times())
    print(location.visibilities(days=3, criterion=1))


if __name__ == "__main__":
    main()
