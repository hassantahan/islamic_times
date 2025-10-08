from datetime import datetime, timedelta
from islamic_times.islamic_times import ITLocation

LOCATIONS = {
    # Australia
    "Brisbane" : (-27.469770, 153.025131, 27),
    "Sydney" : (-33.868820, 151.209290, 58),
    "Canberra" : (-35.280937, 149.130009, 577),
    "Melbourne" : (-37.813629, 144.963058, 31),
    "Adelaide" : (-34.928499, 138.600746, 50),
    "Perth" : (-31.950527, 115.860458, 15),

    # East Asia
    "Tokyo" : (35.689487, 139.691711, 40),
    "Osaka" : (34.693737, 135.502167, 15),
    "Seoul" : (37.566536, 126.977966, 38),
    "Taipei" : (25.033964, 121.564468, 10),
    "Beijing" : (39.904202, 116.407394, 43),
    "Hong Kong" : (22.319303, 114.169361, 9),
    "Singapore" : (1.352083, 103.819839, 15),

    # South Asia
    "Dhaka" : (23.810331, 90.412521, 4),
    "Lucknow" : (26.846695, 80.946167, 123),
    "New Delhi" : (28.613939, 77.209023, 216),
    "Mumbai" : (19.075983, 72.877655, 14),
    "Lahore" : (31.549721, 74.343613, 217),
    "Kabul" : (34.555349, 69.207489, 1790),
    "Hyderabad" : (17.385044, 78.486671, 505),
    "Karachi" : (24.860735, 67.001137, 8),

    # Iran
    "Mashhad" : (36.260462, 59.616755, 995),
    "Kerman" : (30.283937, 57.083362, 1755),
    "Bandar Abbas" : (27.183221, 56.266645, 9),
    "Shiraz" : (29.591768, 52.583698, 1486),
    "Isfahan" : (32.654627, 51.667983, 1575),
    "Qom" : (34.641582, 50.874603, 930),
    "Tehran" : (35.689198, 51.388974, 1191),
    "Tabriz" : (38.096240, 46.273801, 1351),

    # Middle East
    "Muscat" : (23.588000, 58.382900, 15),
    "Dubai" : (25.276987, 55.296249, 5),
    "Abu Dhabi" : (24.453884, 54.377343, 27),
    "Doha" : (25.276987, 51.520008, 10),
    "Bahrain" : (26.066700, 50.557700, 2),
    "Kuwait City" : (29.375859, 47.977405, 15),
    "Basra" : (30.508102, 47.783489, 5),
    "Riyadh" : (24.713552, 46.675297, 612),
    "Najaf" : (31.999820, 44.314730, 60),
    "Karbala" : (32.616550, 44.024996, 30),
    "Baghdad" : (33.315242, 44.366066, 34),
    "Madinah" : (24.524654, 39.569184, 625),
    "Makkah" : (21.389082, 39.857910, 277),
    "Sanaa" : (15.369445, 44.191006, 2200),
    "Erbil" : (36.191113, 44.009167, 420),
    "Mosul" : (36.345581, 43.145279, 223),
    "Aleppo" : (36.202105, 37.134260, 379),
    "Homs" : (34.730837, 36.709408, 501),
    "Damascus" : (33.513807, 36.276528, 680),
    "Amman" : (31.953949, 35.910635, 773),
    "Tripoli, Lebanon" : (34.436668, 35.849724, 10),
    "Beirut" : (33.893791, 35.501778, 50),
    "Sidon" : (33.560635, 35.375740, 10),
    "Tyre" : (33.273709, 35.203766, 5),
    "Jerusalem" : (31.768319, 35.213710, 754),
    "Cairo" : (30.044420, 31.235712, 23),
    "Ankara" : (39.933365, 32.859741, 938),
    "Alexandria" : (31.200092, 29.918739, 5),
    "Istanbul" : (41.008240, 28.978359, 40),

    # Africa
    "Mogadishu" : (2.046934, 45.318161, 9),
    "Djibouti" : (11.572077, 43.145645, 14),
    "Addis Ababa" : (9.019188, 38.752529, 2355),
    "Mombasa" : (-4.043477, 39.668206, 50),
    "Dar es Salaam" : (-6.792354, 39.208328, 53),
    "Khartoum" : (15.500654, 32.559898, 382),
    "Cape Town" : (-33.924870, 18.424055, 25),
    "Tripoli" : (32.887209, 13.191338, 81),
    "Tunis" : (36.806389, 10.181667, 4),
    "Abuja" : (9.057850, 7.495080, 840),
    "Lagos" : (6.524379, 3.379206, 41),
    "Algiers" : (36.753769, 3.058756, 25),
    "Fes" : (34.033127, -5.000280, 415),
    "Rabat" : (34.020882, -6.841650, 75),
    "Marrakesh" : (31.629472, -7.981084, 466),
    "Dakar" : (14.716677, -17.467686, 22),

    # Europe
    "Moscow" : (55.755825, 37.617298, 156),
    "Warsaw" : (52.229676, 21.012229, 100),
    "Stockholm" : (59.329323, 18.068581, 28),
    "Vienna" : (48.210033, 16.373819, 170),
    "Berlin" : (52.520008, 13.404954, 34),
    "Rome" : (41.902782, 12.496366, 21),
    "Oslo" : (59.913868, 10.752245, 23),
    "Milan" : (45.464203, 9.189982, 120),
    "Frankfurt" : (50.110924, 8.682127, 112),
    "Bern" : (46.9481, 7.4474, 540),
    "Amsterdam" : (52.367573, 4.904139, 13),
    "Paris" : (48.856613, 2.352222, 35),
    "London" : (51.507222, -0.127500, 11),
    "Birmingham" : (52.486244, -1.890401, 140),
    "Manchester" : (53.480759, -2.242631, 38),
    "Liverpool" : (53.408371, -2.991573, 70),
    "Madrid" : (40.416775, -3.703790, 667),
    "Lisbon" : (38.716892, -9.139362, 2),

    # NA-East
    "Montreal" : (45.508870, -73.554240, 30),
    "Ottawa" : (45.425230, -75.699960, 54),
    "Toronto" : (43.651070, -79.347015, 150),
    "Hamilton" : (43.255490, -79.873380, 111),
    "London, Canada" : (42.988150, -81.246090, 255),
    "Windsor" : (42.323780, -83.040370, 172),
    "New York" : (40.713050, -74.007230, 105),
    "Philadelphia" : (39.951060, -75.165620, 60),
    "Washington" : (38.892060, -77.019910, 15),
    "Richmond" : (37.540760, -77.433930, 56),
    "Miami" : (25.775080, -80.194700, 16),

    # NA-Central
    "Chicago" : (41.883230, -87.632400, 152),
    "Houston" : (29.760800, -95.369510, 12),
    "Dallas" : (32.777980, -96.796210, 152),
    "Austin" : (30.264980, -97.746600, 151),

    # NA-West
    "Vancouver" : (49.282730, -123.120735, 70),
    "Seattle" : (47.606210, -122.332070, 56),
    "Portland" : (45.505106, -122.675026, 15),
    "San Francisco" : (37.774930, -122.419420, 16),
    "Los Angeles" : (34.052235, -118.243683, 71),
    "San Diego" : (32.715740, -117.161087, 20),

    # South America
    "Rio de Janeiro" : (-22.906847, -43.172897, 5),
    "Sao Paulo" : (-23.550520, -46.633308, 760),
    "Buenos Aires" : (-34.603722, -58.381592, 25),
    "Santiago" : (-33.448891, -70.669266, 520),
    "Mexico City" : (19.432608, -99.133209, 2250),
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
    lat, lon, elev = LOCATIONS["Toronto"]
    dt = datetime.now()
    local: ITLocation = ITLocation(latitude=lat, longitude=lon, elevation=elev, date=dt, find_local_tz=False, method='Jafari')
    main(local, False)