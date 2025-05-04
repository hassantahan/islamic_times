from datetime import datetime
from islamic_times.islamic_times import ITLocation

LOCATIONS = {
    "Montreal, QC" : (45.508870, -73.554240, 30),
    "Ottawa, ON" : (45.425230, -75.699960, 54),
    "Toronto, ON" : (43.651070, -79.347015, 150),
    "Hamilton, ON" : (43.255490, -79.873380, 111),
    "London, ON" : (42.988150, -81.246090, 255),
    "Windsor, ON / Detroit, MI" : (42.323780, -83.040370, 172),
    "New York, NY" : (40.713050, -74.007230, 105),
    "Philadelphia, PA" : (39.951060, -75.165620, 60),
    "Washington, D.C." : (38.892060, -77.019910, 15),
    "Richmond, VA" : (37.540760, -77.433930, 56),
    "Chicago, IL" : (41.883230, -87.632400, 152),
    "Miami, FL" : (25.775080, -80.194700, 16),
    "Houston, TX" : (29.760800, -95.369510, 12),
    "Dallas, TX" : (32.777980, -96.796210, 152),
    "Austin, TX" : (30.264980, -97.746600, 151),
}

def update_time(local: ITLocation):
    local.update_time(datetime.now())
    local.calculate_astro()

def main(local: ITLocation, do_update = True):
    # Observer
    print(local.observer())

    # Astro
    if not local.auto_calculate:
        local.calculate_astro()

    # Date & Time
    print(local.dates_times())

    # Prayer Times
    if not local.auto_calculate:
        local.calculate_prayer_times()
        
    print(local.prayer_times())

    # Mecca
    print(local.mecca())

    # Update Time
    if do_update:
        update_time(local)

    # The Sun
    print(local.sun())

    # The Moon
    print(local.moon())

    # Moon Phases
    print("Moon Phases")
    temp = local.moonphases()
    for phase in temp:
        print(f"\t{phase[0]}\t\t{phase[1].strftime('%X %d-%m-%Y')}")

    # New Moon Visibility
    print(local.visibilities(days=3, criterion=1))

if __name__ == '__main__':    
    lat, lon, elev = LOCATIONS["Toronto, ON"]
    dt = datetime.now()
    local: ITLocation = ITLocation(latitude=lat, longitude=lon, elevation=elev, date=dt, find_local_tz=True)
    main(local, True)
    