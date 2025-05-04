# islamic_times

**islamic_times** is a Python package with an integrated C-extension designed to calculate Islamic astronomical parameters and prayer times. It provides functions for converting Gregorian dates to Hijri, computing sun and moon positions, estimating new moon crescent visibilities, determining Qibla direction, and calculating all major Islamic prayer times according to several methodologies and even allowing for customization of them.

## Features

- **Prayer Times Calculation:**  
  Computes Fajr, Sunrise, Ẓuhr, ʿAṣr (Standard & Ḥanafī), Maghrib, ʿIshāʾ, and Islamic midnight (both types).
  
- **Astronomical Calculations:**  
  Determines solar and lunar positions, including right ascension, declination, altitude, and azimuth.
  
- **Hijri Calendar Conversion:**  
  Converts Gregorian dates to Hijri (Islamic) dates.
  
- **Qibla Direction:**  
  Calculates the distance and bearing (cardinal direction) to Mecca.
  
- **Customizable Methods:**  
  Supports several prayer time calculation methods (e.g., Muslim World League, ISNA, Egypt, Makkah, Karachi, Tehran, Jaʿfarī).
  
- **New Moon Crescent Visibility Calculations**
  Computes the visibility of the nearest new moon for the observer according to either Yallop, 1997 or Odeh, 2006.
  
- **Heavy Computations Done in C**
  For blazing fast performace, all the astronomical calculations are done in an integrated C-extension.
  
- **Extensive Documentation:**  
  Includes detailed API documentation and inline comments for clarity.

## Installation

You can install the package from PyPI using pip:

```bash
pip install islamic_times
```

## API Overview

### ITLocation Class

The core of the package is the `ITLocation` class, which encapsulates both the observer’s parameters and the associated calculations. Key aspects include:

- **Initialization:**
  - `latitude`, `longitude`, `elevation`: Geographic position (i.e. specifically, geodetic). Default is the Greenwich Observatory.
  - `temperature`, `pressure`: Environmental details. Default is 10°C and 101.325 kPA.
  - `date`: Datetime information. Default is current UTC.
  - `method`: Standard Prayer calculation method (e.g., 'JAFARI', 'ISNA', etc.).
  - `find_local_tz`: Automatically find the local timezone. Default is `False` (computationally expensive).
  - `auto_calculate`: Determines if astronomical calculations run immediately. Default is `True`.
  - `asr_type`: Defines the method for ʿAṣr calculation (0 for standard, 1 for Ḥanafī). Default is 0.
  - `midnight_type`: Defines the method for Islamic midnight calculation (0 for middle of sunset → sunrise, 1 for middle of sunset → fajr). Default is 0.

- **Methods:**
  - `update_time(new_date: datetime)`: Updates the observer’s date and time.
  - `calculate_astro()`: Computes astronomical parameters (sun & moon positions, Julian date, etc.). Used if `auto_calculate` is disabled.
  - `calculate_prayer_times()`: Determines prayer times based on the current settings. Used if `auto_calculate` is disabled.
  - `set_prayer_method(method: str, asr_type: int)`: Selects and applies one of the predefined prayer calculation methods.
  - `set_custom_prayer_angles(fajr_angle, maghrib_angle, isha_angle)`: Allows custom adjustment of solar angles.
  - `set_asr_type(asr_type: int)`: Switch between standard and Ḥanafī ʿAṣr calculations.
  - `set_midnight_type(midnight_type: int)`: Customize Islamic midnight calculation.
  - `observer()`: Returns the observer’s geographical/environmental parameters.
  - `dates_times()`: Provides all date and time details, including Hijri conversion.
  - `prayer_times()`: Retrieves the computed prayer times.
  - `mecca()`: Provides distance and Qibla direction to Mecca.
  - `sun()`: Returns information about the Sun’s position and properties.
  - `moon()`: Returns information about the Moon’s position, illumination, and phase.
  - `moonphases()`: Provides details on the nearest moon phases.
  - `visibilities(days: int, criterion: int)`: Computes the visibilities of the nearest new moon crescent for the given amount of days and according to a given criterion.

For a full list of methods and usage details, please refer to the inline documentation within each module.

## Usage

Below is a basic example on how to calculate prayer times for a given location and time:

```python
from datetime import datetime
from islamic_times.islamic_times import ITLocation
from islamic_times.dataclasses import PrayerTimes, Visibilities

location: ITLocation = ITLocation(
    latitude=43.651070,    # Toronto latitude
    longitude=-79.347015,  # Toronto longitude
    elevation=10,          # Elevation in meters
    temperature=15,        # Temperature in °C
    pressure=101.325,      # Atmospheric pressure in kPa
    date=datetime(2005, 6, 1, 12, 0, 0, 0),
    method='ISNA',         # Prayer calculation method
	find_local_tz=True	   # Automatically find the local timezone of the observer
)

# Calculate prayer times and diplay them
prayers: PrayerTimes = location.prayer_times()
print(prayers)

# Calculate the visibilities of the nearest new moon crescent and display them
# (By default, it will compute visibilities for three days from conjunction and according to Yallop, 1997)
vis: Visibilities = location.visibilities()
print(vis)
```

Output:
```
Prayer Times at Observer Timezone
        Method:                 Islamic Society of North America (ISNA)
        Fajr:                   03:52:38 01-06-2005
        Sunrise:                05:38:38 01-06-2005
        Ẓuhr:                   13:15:06 01-06-2005
        ʿAṣr:                   17:19:57 01-06-2005
        Sunset:                 20:52:17 01-06-2005
        Maghrib:                20:52:17 01-06-2005
        ʿIshāʾ:                 22:38:41 01-06-2005
        Midnight:               01:15:13 02-06-2005
Visibility of New Moon Crescent:
        Criterion:              Yallop
        21:05:37 06-06-2005:    -0.753  F: Not visible; below the Danjon limit.
        21:32:17 07-06-2005:    +0.309  A: Easily visible.
        21:54:55 08-06-2005:    +1.517  A: Easily visible.
```

Also included is `test.py` which provides a simple example to showcase the package's functionality.

## Modules

The package is organized into several modules:

- **`islamic_times.py`**:  
  Contains the `ITLocation` class and functions to integrate astronomical calculations with Islamic timings.

- **`it_dataclasses.py`**:  
  Defines core data structures (e.g., `Angle`, `RightAscension`, `Distance`, `ObserverInfo`, `IslamicDateInfo`, and `PrayerMethod`) used across the package.

- **`moon_equations.py`**:  
  Implements lunar calculations including the Moon’s position, phase, and new moon visibility.

- **`prayer_times.py`**:  
  Implements the logic for calculating prayer times using different methodologies and handling extreme latitude cases.

- **`sun_equations.py`**:  
  Contains routines for solar position calculations (declination, right ascension, altitude, azimuth) and the equation of time.

- **`time_equations.py`**:  
  Provides functions for date conversion (Gregorian ⇄ Julian & Hijri), time fractions, sidereal time, and related time equations.

- **`calculation_equations.py`**:  
  Implements various mathematical routines for angle normalization, trigonometric functions in degrees, and other custom calculations.

## Dependencies

The package requires the following external libraries:

- **numpy** – for numerical operations.
- **pytz** – for timezone management.
- **timezonefinder** – (optional) to auto-determine local timezone based on coordinates.
- **dataclasses** – used if running on Python versions earlier than 3.7.

## Default Prayer Calculation Methods

The library supports the following predefined methods:
- **ISNA**: Islamic Society of North America
- **MWL**: Muslim World League
- **Umm al-Qura**: Umm al-Qura University, Makkah
- **Egyptian**: Egyptian General Authority of Survey
- **Karachi**: University of Islamic Sciences, Karachi
- **Tehran**: Institute of Geophysics, University of Tehran 
- **Jafari**: Shia Ithna Ashari, Leva Research Institute, Qom

You can also define custom solar angles for Fajr, Maghrib, and ʿIshāʾ prayers.

---

## Project Structure

```
islamic_times/
├── include/                # Header files for C extensions
│   ├── c_moon_equations.h
│   ├── c_sun_equations.h
│   └── ...
├── src/                    # Source files for C extensions
│   ├── c_moon_equations.c
│   ├── c_sun_equations.c
│   └── ...
├── islamic_times/          # Python modules
│   ├── islamic_times.py    # Main library file
│   ├── prayer_times.py     # Prayer time calculations
│   ├── sun_equations.py    # Sun-related calculations
│   ├── moon_equations.py   # Moon-related calculations
│   └── ...
├── mapper.py               # Mapping script for new moon crescent visibilities
├── test.py                 # Example/testing script
├── setup.py                # Build and installation script
└── README.md               # Project documentation
```

## Mapping New Moon Crescent Visibilities

Included in this package is an efficient mapping tool to map the new moon crescent visibilities for the observer. 
The mapping is done using the `mapper.py` script, which takes the visibility data and generates a visual representation of the crescent visibility across different locations.

The more readily available parameters are found in the `if __name__ == '__main__':` of the script and are as follows:

```python
if __name__ == "__main__":
    # Vars
    today = datetime(2025, 4, 11)
    master_path: str = ""
    total_months: int = 1
    map_region: str = "WORLD" # 'NORTH_AMERICA' 'EUROPE' 'MIDDLE_EAST' 'IRAN' 'WORLD'
    map_mode: str = "category" # "raw" "category"
    resolution: int = 300
    days_to_generate: int = 3
    criterion: int = 1 # Either 0 (Odeh, 2006), or 1 (Yallop, 1997)

    map_region = map_region.upper()
    main()
```

Below is an example of a generated visibility map for the new moon crescent on **2025-04-27 (Dhū al-Qaʿdah 1446)** using the Yallop criterion:

![2025-04-27 Dhū al-Qaʿdah 1446—Yallop](https://github.com/hassantahan/islamic_times/blob/master/2025-04-27%20Dh%C5%AB%20al-Qa%CA%BFdah%201446%E2%80%94Yallop.jpg?raw=true)

## Contributing

All contributions are welcome!

## License

This project is licensed under the [MIT License](LICENSE).

## References

- **Jean Meeus**, *Astronomical Algorithms*, 2nd Edition, Willmann-Bell, Inc., 1998.
- Prayer times calculation methods: [PrayTimes Wiki](http://praytimes.org/wiki/Calculation_Methods).
- Delta T approximation: [NASA's Polynomial Expressions for Delta T (ΔT)](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html)
