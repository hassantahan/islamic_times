# islamic_times: An Astronomical Library for Islamic Astronomical Calculations

The IslamicTimes is a Python library dedicated to providing comprehensive solutions for various Islamic astronomical calculations. This includes but is not limited to, the calculation of prayer times, the positions of the sun and the moon, the moon phases, as well as related astronomical phenomena such as solar and lunar eclipses.

## Features

- **Prayer Times**: The library calculates the five daily Islamic prayer times based on geographical location and date. It includes Fajr, Zuhr, Asr, Maghrib, and Isha.

- **Different Prayer Time Methods Included**: Seven methods to calculate prayer times, along with methods that allow for complete customization of what solar hour angle determines the time of a prayer.

- **Date Conversion**: Conversion between Gregorian and Hijri (Islamic) dates.

- **Sun & Moon Positions**: Provides the position (declination, right ascension, altitude, azimuth) of the sun and the moon at a given time.

- **Moon Illumination**: Calculates the moon's illumination percentage.

- **Moon Phases**: Provides the dates and times of the upcoming moon phases: New Moon, First Quarter, Full Moon, and Last Quarter.

- **New Moon Visibility**: Calculates the visibility of an upcoming new moon according to Yallop, 1997 (a.k.a HMNAO TN No. 69) or Odeh, 2006.

- **Global Lunar Crescent Visibility Maps**: Included is a `mapper.py` script that shows the visibility of a new moon for the day of conjunction, and one and two days after it.

- **Direction to Mecca (Qibla)**: Provides the direction and distance to Mecca from any location worldwide.

## Usage

Below is a simple usage example:

```python
from islamic_times import islamic_times as it
from datetime import datetime

##### Calculation #####
local = it.ITLocation(43.74506, -79.30947, 170.5, datetime.now())

##### Outputs #####
# Date & Time
temp = local.dates_times()
print("Time & Date")
print(f"\tGregorian Date:\t\t{temp['gregorian']}")
print(f"\tIslamic Date:\t\t{temp['hijri']}")
print(f"\t24h-Time:\t\t{temp['time']}\n\tTime Zone:\t\t{temp['timezone']} {temp['utc_offset']}")
print(f"\tLocal JD:\t\t{temp['jd']}")
print(f"\tEquation of time:\t{temp['eq_of_time']} minutes")
print(f"\tEstimated ΔT:\t\t{temp['deltaT']}s")

# Prayer Times
temp = local.calculate_prayer_times()
temp = local.prayer_times()
print("Local Prayer Times")
print(f"\tFajr:\t\t\t{temp['fajr']}")
print(f"\tSunrise:\t\t{temp['sunrise']}")

# More Outputs as required...
```

In the above script, the `ITLocation` class is initialized with the latitude, longitude, and elevation of the location, as well as the date for which calculations are to be performed. The `dates_times()`, `calculate_prayer_times()` and `prayer_times()` methods are then called to retrieve the observer date & time, then calculate the local prayer times (i.e. local according to the observer specified timezone), then retrieve them.

A more comprehensive example of usage can be seen in `test.py`.

## Installation

You can install this package using pip:

```bash
pip install islamic_times
```

## Requirements

Python 3.6 and above is required.

## Contributing

We welcome contributions to this project. Please feel free to open an issue or submit a pull request.

## License

This project is licensed under the MIT License.
