import numpy as np
from time import time
from datetime import datetime, timedelta, timezone
from islamic_times.islamic_times import ITLocation

def decimal_to_dms(decimal_value: float) -> str:
    degrees = int(decimal_value)
    decimal_minutes = (decimal_value - degrees) * 60
    minutes = int(decimal_minutes)
    seconds = (decimal_minutes - minutes) * 60

    dms_string = f"{degrees}° {np.abs(minutes)}' {np.abs(seconds):.2f}\""
    return dms_string

def timing(str_to_print, start_time):
    print(f"{str_to_print}: {(time() - start_time)*1000000:.2f}μs")
    return time()

def output_observer(local: ITLocation):
    start_time = time()
    temp = local.observer()
    timing("Time to fetch observer()", start_time)

    print("Observer Parameters")
    for key, value in temp.items():
        print(f"\t{key.title()}:\t\t{value}")

def calculate_astro(local: ITLocation):
    start_time = time()
    local.calculate_astro()
    timing("Time to calculate_astro()", start_time)

def output_dates_times(local: ITLocation):
    start_time = time()
    temp = local.dates_times()
    timing("Time to fetch dates_times()", start_time)

    print("Time & Date")
    print(f"\tGregorian Date:\t\t{temp['gregorian']}")
    print(f"\tIslamic Date:\t\t{temp['hijri']}")
    print(f"\t24h-Time:\t\t{temp['time']}\n\tTime Zone:\t\t{temp['timezone']} {temp['utc_offset']}")
    print(f"\tLocal JD:\t\t{temp['jd']}")
    print(f"\tEquation of Time:\t{temp['eq_of_time']} minutes")
    print(f"\tEstimated ΔT:\t\t{temp['deltaT']}s")

def update_prayer_times(local: ITLocation):
    # Calculate Prayer Times
    start_time = time()
    local.calculate_prayer_times()
    timing("Time to calculate_prayer_times()", start_time)

def prayer_times(local: ITLocation):
    # Fetch Prayer Times
    start_time = time()
    temp = local.prayer_times()
    timing("Time to fetch prayer_times()", start_time)

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

    # Change Prayer Times Method
    start_time = time()
    local.set_prayer_method('ISNA')
    local.calculate_prayer_times()
    temp = local.prayer_times()
    timing("Time to set new method, recalculate & fetch prayer_times()", start_time)

    print("Changing Prayer Times Method to ISNA")
    print(f"\tMethod:\t\t\t{temp['method']}")
    print(f"\tFajr:\t\t\t{temp['fajr']}")
    print(f"\tMaghrib:\t\t{temp['maghrib']}")
    print(f"\tʿIshāʾ:\t\t\t{temp['isha']}")
    print(f"\tMidnight:\t\t{temp['midnight']}")

    # Custom Prayer Angles
    start_time = time()
    local.set_custom_prayer_angles(fajr_angle=20, maghrib_angle=8, isha_angle=20)
    local.calculate_prayer_times()
    temp = local.prayer_times()
    timing("Time to set custom angles, recalculate & fetch prayer_times()", start_time)

    print("Custom Prayer Times Angles")
    print(f"\tFajr:\t\t\t{temp['fajr']}")
    print(f"\tMaghrib:\t\t{temp['maghrib']}")
    print(f"\tʿIshāʾ:\t\t\t{temp['isha']}")

    # Change ʿAṣr Calculation Method
    start_time = time()
    local.set_asr_type(1)
    local.calculate_prayer_times()
    temp = local.prayer_times()
    timing("Time to change ʿaṣr type, recalculate & fetch prayer_times()", start_time)

    print("Changing ʿAṣr Calculation Method to Ḥanafī")
    print(f"\tʿAṣr:\t\t\t{temp['asr']}\n")

def mecca(local: ITLocation):
    start_time = time()
    temp = local.mecca()
    timing("Time to fetch mecca()", start_time)

    print("Mecca")
    print(f"\tDistance:\t\t{temp['distance']} km")
    print(f"\tDirection:\t\t{temp['cardinal']} ({temp['angle']}°)")

def update_time(local: ITLocation):
    start_time = time()
    local.update_time(datetime.now())
    timing("Time to update_time()", start_time)
    local.calculate_astro()

def output_sun(local: ITLocation):
    start_time = time()
    temp = local.sun()
    timing("Time to fetch sun()", start_time)

    print("The Sun")
    print(f"\tApp. Declination:\t{temp['declination']}°\t{decimal_to_dms(temp['declination'])}")
    print(f"\tApp. Right Ascenscion:\t{temp['right_ascension']}")
    print(f"\tAltitude:\t\t{temp['altitude']}°\t\t{decimal_to_dms(temp['altitude'])}")
    print(f"\tAzimuth:\t\t{temp['azimuth']}°\t\t{decimal_to_dms(temp['azimuth'])}")

def output_moon(local: ITLocation):
    start_time = time()
    temp = local.moon()
    timing("Time to fetch moon()", start_time)

    print("The Moon")
    print(f"\tMoonset:\t\t{temp['moonset']}")
    print(f"\tDeclination:\t\t{temp['declination']}°\t{decimal_to_dms(temp['declination'])}")
    print(f"\tRight Ascenscion:\t{temp['right_ascension']}")
    print(f"\tAltitude:\t\t{temp['altitude']}°\t\t{decimal_to_dms(temp['altitude'])}")
    print(f"\tAzimuth:\t\t{temp['azimuth']}°\t\t{decimal_to_dms(temp['azimuth'])}")
    print(f"\tIllumination:\t\t{temp['illumination']}%")

def output_moonphases(local: ITLocation):
    start_time = time()
    temp = local.moonphases()
    timing("Time to calculate moonphases and fetch moonphases()", start_time)

    print("Moon Phases at Observer Timezone")
    for item in temp:
        print(f"\t{item['phase']}:\t\t{item['datetime'].strftime('%X %A, %d %B, %Y')}")

def output_visibilities(local: ITLocation, days: int = 3, criterion: int = 0):
    start_time = time()
    temp = local.visibilities(days=days, criterion=criterion)
    timing("Time to calculate new moon visibilities and fetch visibilities()", start_time)

    print("Visibility Values of New Moon at Observer TZ 'Best Time' on...")
    for date_label, values in temp.items():
        print(f"\t{date_label}:\t{values[0]:.3f}\t({values[1]})")

def main(local: ITLocation):
    # Observer
    output_observer(local)

    # Astro
    if not local.auto_calculate:
        calculate_astro(local)

    # Date & Time
    output_dates_times(local)

    # Prayer Times
    if not local.auto_calculate:
        update_prayer_times(local)
    prayer_times(local)

    # Mecca
    mecca(local)

    # Update Time
    update_time(local)

    # The Sun
    output_sun(local)

    # The Moon
    output_moon(local)

    # Moon Phases
    output_moonphases(local)

    # New Moon Visibility
    output_visibilities(local)


if __name__ == '__main__':
    start_time = time()
    local = ITLocation()
    timing("Time to initialize", start_time)
    main(local)