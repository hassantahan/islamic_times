from datetime import datetime, timedelta
import numpy as np
from islamic_times import islamic_times as it

def decimal_to_dms(decimal_string: str) -> str:
    # Remove the degree symbol (°) if it's present
    decimal_string = decimal_string.replace("°", "")

    # Convert the input string to a float
    decimal_value = float(decimal_string)

    degrees = int(decimal_value)
    decimal_minutes = (decimal_value - degrees) * 60
    minutes = int(decimal_minutes)
    seconds = (decimal_minutes - minutes) * 60

    dms_string = f"{degrees}° {np.abs(minutes)}' {np.abs(seconds):.2f}\""
    return dms_string

##### Calculation #####
local = it.ITLocation(
                    latitude=34.05224, 
                    longitude=-118.24334, 
                    elevation=120, 
                    temperature=10,
                    pressure=101.325,
                    today=datetime(2024, 12, 28), 
                    find_local_tz=False
                )

##### Outputs #####
# Observer
temp = local.observer()
print("Observer Parameters")
for key, value in temp.items():
    print(f"\t{key}:\t\t{value}")

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
local.calculate_prayer_times()
temp = local.prayer_times()
print("Prayer Times at Observer Timezone")
print(f"\tMethod:\t\t\t{temp['method']}")
print(f"\tFajr:\t\t\t{temp['fajr']}")
print(f"\tSunrise:\t\t{temp['sunrise']}")
print(f"\tẒuhr:\t\t\t {temp['noon']}")
print(f"\tʿAṣr:\t\t\t{temp['asr']}")
print(f"\tSunset:\t\t\t{temp['sunset']}")
print(f"\tMaghrib:\t\t{temp['maghrib']}")
print(f"\tʿIshāʾ:\t\t\t{temp['isha']}")
print(f"\tMidnight:\t\t{temp['midnight']}")

# Mecca
temp = local.mecca()
print("Mecca")
print(f"\tDistance:\t\t{temp['distance']} km")
print(f"\tDirection:\t\t{temp['cardinal']} ({temp['angle']}°)")

# The Sun
temp = local.sun()
print("The Sun")
print(f"\tApp. Declination:\t{temp['declination']}\t{decimal_to_dms(temp['declination'])}")
print(f"\tApp. Right Ascenscion:\t{temp['right_ascension']}")
print(f"\tAltitude:\t\t{temp['altitude']}\t\t{decimal_to_dms(temp['altitude'])}")
print(f"\tAzimuth:\t\t{temp['azimuth']}\t\t{decimal_to_dms(temp['azimuth'])}")

# The Moon
temp = local.moon()
print("The Moon")
print(f"\tMoonset:\t\t{temp['moonset']}")
print(f"\tDeclination:\t\t{temp['declination']}\t{decimal_to_dms(temp['declination'])}")
print(f"\tRight Ascenscion:\t{temp['right_ascension']}")
print(f"\tAltitude:\t\t{temp['altitude']}\t\t{decimal_to_dms(temp['altitude'])}")
print(f"\tAzimuth:\t\t{temp['azimuth']}\t\t{decimal_to_dms(temp['azimuth'])}")
print(f"\tIllumination:\t\t{temp['illumination']}")

# Moon Phases
temp = local.moonphases()
print("Moon Phases at Observer Timezone")
for item in temp:
    print(f"\t{item['phase']}:\t\t{item['datetime'].strftime('%H:%M:%S %A, %d %B, %Y')}")

# New Moon Visibility
temp = local.visibilities(type=1)
print("Visibility Values of New Moon at Observer TZ 'Best Time' on...")
for date_label, values in temp.items():
    print(f"\t{date_label}:\t{values[0]:.3f}\t({values[1]})")