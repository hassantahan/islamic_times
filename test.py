from islamic_times import islamic_times as it
from datetime import datetime, timedelta

##### Definitions #####
TO_LAT 	= 43.74506
TO_LONG = -79.30947
TO_ELEV = 170.5

##### Inputs #####
latitude = TO_LAT #37.336111
longitude = TO_LONG #-121.890556
elev = TO_ELEV #25

##### Calculation #####
local = it.ITLocation(TO_LAT, TO_LONG, TO_ELEV, datetime.now() + timedelta(days=20))

##### Outputs #####
# Date & Time
print("Time & Date\n\tGregorian Date:\t\t{}".format(local.datetime()["gregorian"]))
print(f"\tIslamic Date:\t\t{local.datetime()['hijri']}")
print("\t24h-Time:\t\t{}\n\tTime Zone:\t\t{} {}".format(local.datetime()["time"], local.datetime()["timezone"], local.datetime()["utc_offset"]))
print("\tEquation of time:\t{} minutes".format(local.datetime()["eq_of_time"]))

# Prayer Times
print("Prayer Times\n\tFajr:\t\t\t{}".format(local.prayertimes()["fajr"]))
print("\tSunrise:\t\t{}".format(local.prayertimes()["sunrise"]))
print("\tẒuhr:\t\t\t {}".format(local.prayertimes()["noon"]))
print("\t‘Asr:\t\t\t{}".format(local.prayertimes()["asr"]))
print("\tSunset: \t\t{}".format(local.prayertimes()["sunset"]))
print("\tMaghrib: \t\t{}".format(local.prayertimes()["maghrib"]))
print("\t‘Isha: \t\t\t{}".format(local.prayertimes()["isha"]))
print("\tMidnight: \t\t{}".format(local.prayertimes()["midnight"]))

# Mecca
print("Mecca\n\tDistance: \t\t{} km".format(local.mecca()["distance"]))
print("\tDirection: \t\t{} ({}°)".format(local.mecca()["cardinal"], local.mecca()["angle"]))

# The Sun
print("The Sun\n\tApp. Declination:\t{}".format(local.sun()["declination"]))
print("\tApp. Right Ascenscion:\t{}".format(local.sun()["right_ascension"]))
print("\tAltitude:\t\t{}".format(local.sun()["altitude"]))
print("\tAzimuth:\t\t{}".format(local.sun()["azimuth"]))

# The Moon
print("The Moon\n\tApp. Declination:\t{}".format(local.moon()["declination"]))
print("\tApp. Right Ascenscion:\t{}".format(local.moon()["right_ascension"]))
print("\tAltitude:\t\t{}".format(local.moon()["altitude"]))
print("\tAzimuth:\t\t{}".format(local.moon()["azimuth"]))
print("\tIllumination:\t\t{}".format(local.moon()["illumination"]))

# Moon Phases
print("Moon Phases\n\t{}:\t\t{}".format(local.moonphases()[0]["phase"], local.moonphases()[0]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(local.moonphases()[1]["phase"], local.moonphases()[1]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(local.moonphases()[2]["phase"], local.moonphases()[2]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(local.moonphases()[3]["phase"], local.moonphases()[3]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))

# New Moon Visibility
print("Visibility Values of New Moon after\n\t0 days:\t\t\t{:.3f} ({})".format(local.visibilities()["0"][0], local.visibilities()["0"][1]))
print("\t1 day:\t\t\t{:.3f} ({})".format(local.visibilities()["1"][0], local.visibilities()["1"][1]))
print("\t2 days:\t\t\t{:.3f} ({})".format(local.visibilities()["2"][0], local.visibilities()["2"][1]))