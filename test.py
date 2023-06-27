import islamic_times
from datetime import datetime

##### Definitions #####
TO_LAT 	= 43.74506
TO_LONG = -79.30947
TO_ELEV = 170.5

##### Inputs #####
latitude = TO_LAT #37.336111
longitude = TO_LONG #-121.890556
elev = TO_ELEV #25

##### Calculation #####
it = islamic_times.ITLocation(TO_LAT, TO_LONG, TO_ELEV, datetime(2076, 11, 28, 11, 11, 11))

##### Outputs #####
# Date & Time
print("Time & Date\n\tGregorian Date:\t\t{}".format(it.datetime()["gregorian"]))
print(f"\tIslamic Date:\t\t{it.datetime()['hijri']}")
print("\t24h-Time:\t\t{}\n\tTime Zone:\t\t{} {}".format(it.datetime()["time"], it.datetime()["timezone"], it.datetime()["utc_offset"]))
print("\tEquation of time:\t{} minutes".format(it.datetime()["eq_of_time"]))

# Prayer Times
print("Prayer Times\n\tFajr:\t\t\t{}".format(it.prayertimes()["fajr"]))
print("\tSunrise:\t\t{}".format(it.prayertimes()["sunrise"]))
print("\tẒuhr:\t\t\t {}".format(it.prayertimes()["noon"]))
print("\t‘Asr:\t\t\t{}".format(it.prayertimes()["asr"]))
print("\tSunset: \t\t{}".format(it.prayertimes()["sunset"]))
print("\tMaghrib: \t\t{}".format(it.prayertimes()["maghrib"]))
print("\t‘Isha: \t\t\t{}".format(it.prayertimes()["isha"]))
print("\tMidnight: \t\t{}".format(it.prayertimes()["midnight"]))

# Mecca
print("Mecca\n\tDistance: \t\t{} km".format(it.mecca()["distance"]))
print("\tDirection: \t\t{} ({}°)".format(it.mecca()["cardinal"], it.mecca()["angle"]))

# The Sun
print("The Sun\n\tApp. Declination:\t{}".format(it.sun()["declination"]))
print("\tApp. Right Ascenscion:\t{}".format(it.sun()["right_ascension"]))
print("\tAltitude:\t\t{}".format(it.sun()["altitude"]))
print("\tAzimuth:\t\t{}".format(it.sun()["azimuth"]))

# The Moon
print("The Moon\n\tApp. Declination:\t{}".format(it.moon()["declination"]))
print("\tApp. Right Ascenscion:\t{}".format(it.moon()["right_ascension"]))
print("\tAltitude:\t\t{}".format(it.moon()["altitude"]))
print("\tAzimuth:\t\t{}".format(it.moon()["azimuth"]))
print("\tIllumination:\t\t{}".format(it.moon()["illumination"]))

# Moon Phases
print("Moon Phases\n\t{}:\t\t{}".format(it.moonphases()[0]["phase"], it.moonphases()[0]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(it.moonphases()[1]["phase"], it.moonphases()[1]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(it.moonphases()[2]["phase"], it.moonphases()[2]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(it.moonphases()[3]["phase"], it.moonphases()[3]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))

# New Moon Visibility
print("Visibility Values of New Moon after\n\t0 days:\t\t\t{:.3f} ({})".format(it.visibilities()["0"][0], it.visibilities()["0"][1]))
print("\t1 day:\t\t\t{:.3f} ({})".format(it.visibilities()["1"][0], it.visibilities()["1"][1]))
print("\t2 days:\t\t\t{:.3f} ({})".format(it.visibilities()["2"][0], it.visibilities()["2"][1]))