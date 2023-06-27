# islamic_times: An Astronomical Library for Islamic Astronomical Calculations

The IslamicTimes is a Python library dedicated to providing comprehensive solutions for various Islamic astronomical calculations. This includes but is not limited to, the calculation of prayer times, the positions of the sun and the moon, the moon phases, as well as related astronomical phenomena such as solar and lunar eclipses.

## Features

- **Prayer Times**: The library calculates the five daily Islamic prayer times based on geographical location and date. It includes Fajr, Zuhr, Asr, Maghrib, and Isha.

- **Date Conversion**: Conversion between Gregorian and Hijri (Islamic) dates.

- **Sun & Moon Positions**: Provides the position (declination, right ascension, altitude, azimuth) of the sun and the moon at a given time.

- **Moon Illumination**: Calculates the moon's illumination percentage.

- **Moon Phases**: Provides the dates and times of the upcoming moon phases: New Moon, First Quarter, Full Moon, and Last Quarter.

- **New Moon Visibility**: Calculates the visibility of the new moon, which is crucial for determining the beginning of Islamic months.

- **Direction to Mecca (Qibla)**: Provides the direction and distance to Mecca from any location worldwide.

## Usage

Below is a simple usage example:

```python
import islamic_times
from datetime import datetime

# Definitions
TO_LAT 	= 43.745
TO_LONG = -79.309
TO_ELEV = 150

# Calculation
it = islamic_times.ITLocation(TO_LAT, TO_LONG, TO_ELEV, datetime(2023, 6, 27, 11, 10, 52))

# Outputs
print("Time & Date\n\tGregorian Date:\t\t{}".format(it.datetime()["gregorian"]))
print(f"\tIslamic Date:\t\t{it.datetime()['hijri']}")

# Prayer Times
print("Prayer Times\n\tFajr:\t\t\t{}".format(it.prayertimes()["fajr"]))
print("\tẒuhr:\t\t\t {}".format(it.prayertimes()["noon"]))

# More Outputs as required...
```

In the above script, the `ITLocation` class is initialized with the latitude, longitude, and elevation of the location, as well as the date for which calculations are to be performed. The `datetime` and `prayertimes` methods are then called to retrieve the date and prayer times respectively.

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
