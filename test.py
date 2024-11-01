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
latitude = TO_LAT #37.336111
longitude = TO_LONG #-121.890556
elev = TO_ELEV #25

##### Calculation #####
local = it.ITLocation(TO_LAT, TO_LONG, TO_ELEV, datetime.now())#datetime(2024, 3, 1) )#, datetime(year=2002, month=2, day=12, hour=16, minute=19, second=12))

##### Outputs #####
# Date & Time
temp = local.datetime()
print("Time & Date\n\tGregorian Date:\t\t{}".format(temp["gregorian"]))
print(f"\tIslamic Date:\t\t{temp['hijri']}")
print("\t24h-Time:\t\t{}\n\tTime Zone:\t\t{} {}".format(temp["time"], temp["timezone"], temp["utc_offset"]))
print("\tEquation of time:\t{} minutes".format(temp["eq_of_time"]))
print("\tEstimated ΔT:\t\t{}s".format(temp["deltaT"]))

# Prayer Times
temp = local.prayertimes()
print("Local Prayer Times\n\tFajr:\t\t\t{}".format(temp["fajr"]))
print("\tSunrise:\t\t{}".format(temp["sunrise"]))
print("\tẒuhr:\t\t\t {}".format(temp["noon"]))
print("\tʿAsr:\t\t\t{}".format(temp["asr"]))
print("\tSunset: \t\t{}".format(temp["sunset"]))
print("\tMaghrib: \t\t{}".format(temp["maghrib"]))
print("\tʿIsha: \t\t\t{}".format(temp["isha"]))
print("\tMidnight: \t\t{}".format(temp["midnight"]))

# Mecca
temp = local.mecca()
print("Mecca\n\tDistance: \t\t{} km".format(temp["distance"]))
print("\tDirection: \t\t{} ({}°)".format(temp["cardinal"], temp["angle"]))

# The Sun
temp = local.sun()
print("The Sun\n\tApp. Declination:\t{}\t\t{}".format(temp["declination"], decimal_to_dms(temp["declination"])))
print("\tApp. Right Ascenscion:\t{}".format(temp["right_ascension"]))
print("\tAltitude:\t\t{}\t\t{}".format(temp["altitude"], decimal_to_dms(temp["altitude"])))
print("\tAzimuth:\t\t{}\t\t{}".format(temp["azimuth"], decimal_to_dms(temp["azimuth"])))

# The Moon
temp = local.moon()
print("The Moon\n\tApp. Declination:\t{}\t{}".format(temp["declination"], decimal_to_dms(temp["declination"])))
print("\tApp. Right Ascenscion:\t{}".format(temp["right_ascension"]))
print("\tAltitude:\t\t{}\t\t{}".format(temp["altitude"], decimal_to_dms(temp["altitude"])))
print("\tAzimuth:\t\t{}\t\t{}".format(temp["azimuth"], decimal_to_dms(temp["azimuth"])))
print("\tIllumination:\t\t{}".format(temp["illumination"]))

# Moon Phases
temp = local.moonphases()
print("Moon Phases\n\t{}:\t\t{}".format(temp[0]["phase"], temp[0]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(temp[1]["phase"], temp[1]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(temp[2]["phase"], temp[2]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(temp[3]["phase"], temp[3]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))

# New Moon Visibility
temp = local.visibilities()
print("Visibility Values of New Moon after...\n\t0 days:\t\t\t{:.3f} ({})".format(temp["0"][0], temp["0"][1]))
print("\t1 day:\t\t\t{:.3f} ({})".format(temp["1"][0], temp["1"][1]))
print("\t2 days:\t\t\t{:.3f} ({})".format(temp["2"][0], temp["2"][1]))