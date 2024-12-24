from islamic_times import islamic_times as it
from datetime import datetime, timedelta
import numpy as np

def decimal_to_dms(decimal_string):
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

##### Definitions #####
TO_LAT 	= 43.74506
TO_LONG = -79.30947
TO_ELEV = 170.5

HAM_LAT = 43.232737
HAM_LONG = -79.857990
HAM_ELEV = 80

MON_LAT = 45.53932
MON_LONG = -73.58525
MON_ELEV = 50

##### Inputs #####
latitude = TO_LAT
longitude = TO_LONG
elev = TO_ELEV

##### Calculation #####
local = it.ITLocation(TO_LAT, TO_LONG, TO_ELEV, datetime.now() + timedelta(days=90))

##### Outputs #####
# Date & Time
temp = local.datetime()
print("Time & Date")
print(f"\tGregorian Date:\t\t{temp['gregorian']}")
print(f"\tIslamic Date:\t\t{temp['hijri']}")
print(f"\t24h-Time:\t\t{temp['time']}\n\tTime Zone:\t\t{temp['timezone']} {temp['utc_offset']}")
print(f"\tLocal JD:\t\t{temp['jd']}")
print(f"\tEquation of time:\t{temp['eq_of_time']} minutes")
print(f"\tEstimated ΔT:\t\t{temp['deltaT']}s")

# Prayer Times
temp = local.prayertimes()
print("Local Prayer Times")
print(f"\tFajr:\t\t\t{temp['fajr']}")
print(f"\tSunrise:\t\t{temp['sunrise']}")
print(f"\tẒuhr:\t\t\t {temp['noon']}")
print(f"\tʿAsr:\t\t\t{temp['asr']}")
print(f"\tSunset:\t\t\t{temp['sunset']}")
print(f"\tMaghrib:\t\t{temp['maghrib']}")
print(f"\tʿIsha:\t\t\t{temp['isha']}")
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
print(f"\tApp. Right Ascenscion:\t{temp['right_ascension']}".format())
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
print("Moon Phases")
print(f"\t{temp[0]['phase']}:\t\t{temp[0]['datetime'].strftime('%H:%M:%S %A, %d %B, %Y')}")
print(f"\t{temp[1]['phase']}:\t\t{temp[1]['datetime'].strftime('%H:%M:%S %A, %d %B, %Y')}")
print(f"\t{temp[2]['phase']}:\t\t{temp[2]['datetime'].strftime('%H:%M:%S %A, %d %B, %Y')}")
print(f"\t{temp[3]['phase']}:\t\t{temp[3]['datetime'].strftime('%H:%M:%S %A, %d %B, %Y')}")

# New Moon Visibility
temp = local.visibilities()
print("Visibility Values of New Moon after...")
print(f"\t0 days:\t\t\t{temp['0'][0]:.3f} ({temp['0'][1]})")
print(f"\t1 day:\t\t\t{temp['1'][0]:.3f} ({temp['1'][1]})")
print(f"\t2 days:\t\t\t{temp['2'][0]:.3f} ({temp['2'][1]})")