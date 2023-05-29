import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib
from main_2 import *
from mpl_toolkits.basemap import Basemap

def calculate_q_value(date, time, lat, long):
    arcv, arcl, daz, pi, alt = calculate_values(date, time, lat, long)
    
    semi_diameter = 0.27245 * pi
    semi_diameter_prime = semi_diameter * (1 + math.sin(math.radians(alt)) * math.sin(math.radians(pi / 60)))

    w_prime = semi_diameter_prime * (1 - math.cos(math.radians(arcl)))

    q_value = (arcv - (11.8371 - 6.3226 * w_prime + 0.7319 * w_prime ** 2 - 0.1018 * w_prime ** 3)) / 10

    if q_value > 10:
        print(f'date: {date}, time: {time}, q_value: {q_value:.2f}, arcv: {arcv:.2f}, arcl: {arcl:.2f}, daz: {daz:.2f} lat: {lat}, long: {long}')
        q_value = -np.inf

    return q_value

def generate_map(calculating_date, day_of_new_moon):
    # Define the latitudes and longitudes with a resolution of 1 degree
    latitudes = np.arange(-55, 66, 1)
    longitudes = np.arange(-179, 181, 1)

    q_values = np.empty((len(latitudes), len(longitudes)))
    #classified_values = np.empty((len(latitudes), len(longitudes)))
    
    for i, lat in enumerate(latitudes):
        for j, long in enumerate(longitudes):
            # Set up local observer
            observer = ephem.Observer()
            observer.lat = str(lat)
            observer.long = str(long)
            
            utc_offset = get_utc_offset(lat, long, calculating_date)

            k = calculating_date - timedelta(hours=utc_offset)
            observer.date = f'{k:%Y/%m/%d %H:%M:%S}'

            sunset_utc = observer.next_setting(ephem.Sun())
            sunset_local = convert_utc_to_offset(str(sunset_utc), utc_offset)

            new_moon_datetime_utc = calculating_date.strftime('%Y/%m/%d %H:%M:%S')
            new_moon_datetime_local = convert_utc_to_offset(new_moon_datetime_utc, utc_offset)

            d1 = datetime.strptime(sunset_local, '%Y/%m/%d %H:%M:%S')
            d2 = datetime.strptime(new_moon_datetime_local, '%Y/%m/%d %H:%M:%S')

            sunset_after_new_moon = 0

            if d1 > d2:
                #Sunset happens AFTER the new moon.
                sunset_after_new_moon = 1

            if sunset_after_new_moon or not(day_of_new_moon):
                sunset_utc_datetime = datetime.strptime(str(sunset_utc), '%Y/%m/%d %H:%M:%S')

                # Format the date and time as strings
                date_str = sunset_utc_datetime.strftime('%Y/%m/%d')
                time_str = sunset_utc_datetime.strftime('%H:%M:%S')
                q_values[i, j] = calculate_q_value(date_str, time_str, lat, long)                                
            else:
                q_values[i,j] = -np.inf
    
    # Generate map
    plt.figure(figsize=(9, 5))

    # Create a new Basemap instance with the desired projection
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)

    # Draw coastlines
    m.drawcoastlines()

    # Plot the pcolormesh on the Basemap instance
    c = m.pcolormesh(longitudes, latitudes, q_values, cmap='viridis')

    # Filter out -np.inf values from q_values
    valid_values = q_values[q_values > -np.inf]

    # Add value lines
    num_lines = 5  # Number of value lines to add
    min_value = np.min(valid_values)  # Minimum value in the q_values array
    max_value = np.max(valid_values)  # Maximum value in the q_values array
    value_step = (max_value - min_value) / (num_lines + 1)  # Step size between value lines
    value_lines = np.linspace(min_value + value_step, max_value - value_step, num_lines)  # Array of value lines

    # Reshape longitudes and latitudes to 2D arrays
    lon_2d, lat_2d = np.meshgrid(longitudes, latitudes)

    # Draw contour lines with labels
    contour = m.contour(lon_2d, lat_2d, q_values, levels=value_lines, colors='red')
    plt.clabel(contour, inline=True, fmt='%.2f', colors='red', fontsize=8)  # Label format and style

    # Add a colorbar with tick labels representing the classified values
    cbar = plt.colorbar(c, label='Crescent Visibility')

    # Add a legend specifying the contour line values
    #legend_labels = ['%.2f' % value for value in value_lines]
    #plt.legend(contour.collections, legend_labels, title='Contour Lines', loc='lower right')

    # Add a title
    k = calculating_date.strftime('%Y/%m/%d')
    plt.title(f'Crescent Visibility on {k}')

    # Show the plot
    plt.show()

    # Define conditions and corresponding values
    conditions = [
        (q_values > 0.216), 
        (q_values <= 0.216) & (q_values > -0.014), 
        (q_values <= -0.014) & (q_values > -0.160), 
        (q_values <= -0.160) & (q_values > -0.232), 
        (q_values <= -0.232) & (q_values > -0.293), 
        (q_values <= -0.293)
    ]
    values = [0, 1, 2, 3, 4, 5]  # Mapping each condition to an integer

    # Apply conditions to q_values array
    q_levels = np.select(conditions, values)

    plt.figure(figsize=(12, 7))
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines()

    # Define colormap
    cmap = matplotlib.colors.ListedColormap(['green', 'blue', 'yellow', 'orange', 'red', 'purple'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    # Plot the pcolormesh on the Basemap instance with the new color array
    c = m.pcolormesh(longitudes, latitudes, q_levels, cmap=cmap, norm=norm)

    # Create a custom legend
    legend_labels = ['Crescent easily visible', 'Crescent visible under perfect conditions', 'May need optical aid to find crescent', 
                     'Will need optical aid to find crescent', 'Crescent not visible with telescope', 'Crescent not visible, below the Danjon limit']
    legend_colors = [matplotlib.patches.Patch(color=clr, label=lbl) for clr, lbl in zip(cmap.colors, legend_labels)]
    plt.legend(handles=legend_colors, loc='lower right')

    plt.title(f'Crescent Visibility on {k}')
    plt.show()

current_date = datetime.now() #datetime(2023, 5, 18, 0, 0, 0, 0, tzinfo=timezone.utc)
next_new_moon_date = str(next_new_moon(current_date))
dt_next_new_moon_date = datetime.strptime(next_new_moon_date, '%Y/%m/%d %H:%M:%S')
print(f'New Moon occurs on: {dt_next_new_moon_date} UTC')
# generate_map(dt_next_new_moon_date + timedelta(days=-1), 1)
generate_map(dt_next_new_moon_date, 1)
generate_map(dt_next_new_moon_date + timedelta(days=1), 0)
generate_map(dt_next_new_moon_date + timedelta(days=2), 0)