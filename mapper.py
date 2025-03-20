import os
import sys
import numpy as np
from time import time
import geopandas as gpd
from datetime import timedelta, datetime
from islamic_times import islamic_times as it
from islamic_times.time_equations import get_islamic_month, gregorian_to_hijri

CITIES = [
        'Washington,  D.C.', 'Los Angeles', 'London',
        'Vienna', 'Moscow', 'Tokyo', 'Beijing', 'Mumbai',
        'Cairo', 'Lagos', 'Sydney', 'São Paulo', 'Mexico City', 
        'Toronto', 'Istanbul', 'Tehran', 'Islamabad', 'Perth',
        'Bangkok', 'Hong Kong', 'Singapore', 'Makkah', 'Miami',
        'Lima', 'Bogota', 'Cape Town', 'Madrid', 'Vancouver',
        'Nairobi', 'Addis Ababa', 'Dakar', 'Santiago'
    ]

class Tee:
    def __init__(self, filename, mode="w", encoding="utf-8"):
        self.file = open(filename, mode, encoding=encoding)
        self.stdout = sys.stdout

    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)

    def flush(self):
        self.stdout.flush()
        self.file.flush()

def compute_q_values_for_lat(i_lat: float, lon: float, day: datetime, amount: int, type: int):
        q_values_temp = []
        for i_lon in lon:
            local = it.ITLocation(latitude=i_lat, longitude=i_lon, elevation=0, date=day, find_local_tz=False, auto_calculate=False)
            vis = local.visibilities(days=amount, criterion=type)
            for label, value in vis.items():
                q_values_temp.append(value[1])
        return q_values_temp

def calculate_moon_visibility(lon: float, lat: float, day: datetime, amount: int = 1, type: int = 0):
    moon_phases = it.ITLocation(date=day, auto_calculate=False).moonphases()
    for item in moon_phases:
            if item['phase'] == "New Moon":
                new_moon_date = item['datetime']
                break  
    
    from multiprocessing import Pool
    def parallel_compute_q_values(lat: float, lon: float, day: datetime, amount: int, type: int, num_processes: int = 8):
        # Create the pool of worker processes
        with Pool(processes=num_processes) as pool:
            # Map each latitude to the compute_q_values_for_lat function
            results = pool.starmap(
                compute_q_values_for_lat,
                [(i_lat, lon, day, amount, type) for i_lat in lat]  # arguments per latitude
            )
        return results

    q_values = parallel_compute_q_values(lat, lon, day, amount, type)

    return q_values, new_moon_date

def load_shapefiles(states_path, places_path, cities):
    states_gdf = gpd.read_file(states_path)
    places_gdf = gpd.read_file(places_path)
    places_gdf = places_gdf[places_gdf['NAME'].isin(cities)]
    places_gdf = places_gdf.loc[places_gdf.groupby('NAME')['POP_MAX'].idxmax()]
    return states_gdf, places_gdf

def clip_map(states_gdf, places_gdf, minx=-179, maxx=180, miny=-61, maxy=61):
    from shapely.geometry import Polygon
    bbox_polygon = Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)])
    bbox_gdf = gpd.GeoDataFrame([1], geometry=[bbox_polygon], crs=states_gdf.crs)
    states_clip = gpd.overlay(states_gdf, bbox_gdf, how='intersection')
    places_clip = places_gdf.cx[minx:maxx, miny:maxy]
    return states_clip, places_clip

def create_grid(resolution, minx=-179, maxx=180, miny=-61, maxy=61):
    lon_vals = np.linspace(minx, maxx, resolution)
    lat_vals = np.linspace(miny, maxy, resolution)
    return lon_vals, lat_vals, len(lon_vals), len(lat_vals)

def get_category_colors(visibility_type):
    import matplotlib.colors as mcolors
    categories = {
        0: {
            "Moonset before the new moon.": "#141414",
            "Moonset before sunset.": "#393a3c",
            "D: Crescent is not visible even by optical aid.": "#807f80",
            "C: Crescent is visible by optical aid only.": "#B89D18",
            "B: Crescent is visible by optical aid, and it could be seen by naked eyes.": "#74b818",
            "A: Crescent is visible by naked eyes.": "#1BB818"
        },
        1: {
            "Moonset before the new moon.": "#141414",
            "Moonset before sunset.": "#393a3c",
            "F: Not visible; below the Danjon limit.": "#807f80",
            "E: Not visible with a [conventional] telescope.": "#B81818",
            "D: Will need optical aid to find crescent.": "#e3d61b",
            "C: May need optical aid to find crescent.": "#89d518",
            "B: Visible under perfect conditions.": "#54b818",
            "A: Easily visible.": "#1bdf18",
        }
    }
    selected_categories = categories[visibility_type]
    transparent_category = list(selected_categories.keys())[2]
    colors_rgba = {k: mcolors.to_rgba(v) if k != transparent_category else (0, 0, 0, 0.1) for k, v in selected_categories.items()}
    return selected_categories, colors_rgba

def map_visibilities(visibilities_3d, category_to_index, ny, nx, amount):
    visibilities_mapped = np.zeros((ny, nx, amount), dtype=int)
    for i_lat in range(ny):
        for i_lon in range(nx):
            for i_day in range(amount):
                visibilities_mapped[i_lat, i_lon, i_day] = category_to_index[visibilities_3d[i_lat, i_lon, i_day]]
    return visibilities_mapped

def plot_map(lon_vals, lat_vals, visibilities_mapped, states_clip, places_clip, unique_categories, category_colors_rgba, start_date, amount, out_dir):
    from textwrap import wrap
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Rectangle
    import matplotlib.patheffects as path_effects
    
    cmap = mcolors.ListedColormap([category_colors_rgba[cat] for cat in unique_categories])
    bounds = np.arange(len(unique_categories) + 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    fig = plt.figure(figsize=(20, 15), dpi=80)
    gs = gridspec.GridSpec(amount, 2, width_ratios=[50, 1], height_ratios=[2] * amount)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(amount)]

    for i_day, ax in enumerate(axes):
        mesh = ax.pcolormesh(lon_vals, lat_vals, visibilities_mapped[..., i_day], cmap=cmap, norm=norm, shading="auto")
        states_clip.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=0.5)
        places_clip.plot(ax=ax, color='violet', markersize=7)
        for idx, row in places_clip.iterrows():
            text = ax.text(row.geometry.x + 0.05, row.geometry.y + 0.05, row['NAME'], fontsize=12, color='white', ha='left', va='bottom')
            text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'), path_effects.Normal()])
        ax.set_xlim(min(lon_vals), max(lon_vals))
        ax.set_ylim(min(lat_vals), max(lat_vals))
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"New Moon Visibility on {(start_date + timedelta(days=i_day)).strftime('%d-%m-%Y')} at Local Best Time")
    
    # **Creating the Legend**
    legend_ax = fig.add_subplot(gs[:, 1])  # Assign the second column for legend
    legend_ax.axis("off")  # Hide axes

    def wrap_text(text, width=20):
        """Wrap text into multiple lines for better display in the legend."""
        return "\n".join(wrap(text, width))

    # Adjust spacing between legend items
    legend_labels = list(unique_categories)
    row_height = 1.05  # Increased row spacing to prevent text overlap

    for idx, category in enumerate(legend_labels):
        wrapped_category = wrap_text(category, width=30)  # Wrap long text
        legend_ax.add_patch(Rectangle((0, idx * row_height), 1, 1, color=category_colors_rgba[category]))
        legend_ax.text(3.2, idx * row_height + 0.5, wrapped_category, fontsize=12, va='center', ha='left')

    legend_ax.set_xlim(0, 2)
    legend_ax.set_ylim(0, len(legend_labels) * row_height)
    legend_ax.set_aspect("auto")

    criterion_string = "Odeh, 2006" if criterion == 0 else "Yallop, 1997"

    # Adjust the plots slightly
    plt.subplots_adjust(hspace=0.2, left=0.05, right=0.85, top=0.95, bottom=0.05)

    # Add text
    plt.figtext(0.15, 0.01, f"The New Moon (i.e. conjunction) occurs at {start_date.strftime("%Y-%m-%d %X")} UTC", ha="center", fontsize=12)
    plt.figtext(0.945, 0.03, f"Criterion: {criterion_string}", ha="center", fontsize=12)
    plt.figtext(0.84, 0.01, f"CC BY-SA | Hassan Tahan | Created with the islamic_times Python library", ha="center", fontsize=12)
    plt.figtext(0.5, 0.98, f"{days}-Day New Moon Crescent Visibility Map for {islamic_month_name}, {islamic_year} A.H.", ha="center", fontsize=16)

    # Name the file
    name: str = f"{start_date.strftime("%Y-%m-%d")} {islamic_month_name} {islamic_year}"
    if criterion == 1:
        name += " Type 1"

    # Save and close
    plt.savefig(os.path.join(out_dir, name))
    plt.close()

def print_ts(message: str):
    print(f"[{datetime.fromtimestamp(time()).strftime("%X %d-%m-%Y")}] {message}")

def main(day, res, amount=1, visibility_criterion=0, coords=(-179, 180, -61, 61), path=""):
    print_ts(f"Loading shape files...")
    t1 = time()
    states_gdf, places_gdf = load_shapefiles("map_shp_files/combined_polygons.shp", "map_shp_files/combined_points.shp", CITIES)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    print_ts(f"Clipping map to coordinates...")
    t1 = time()
    states_clip, places_clip = clip_map(states_gdf, places_gdf, minx=coords[0], maxx=coords[1], miny=coords[2], maxy=coords[3])
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    print_ts(f"Creating map grid...")
    t1 = time()
    lon_vals, lat_vals, nx, ny = create_grid(res, minx=coords[0], maxx=coords[1], miny=coords[2], maxy=coords[3])
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    print_ts(f"Calculating new moon crescent visibilities...")
    t1 = time()
    visibilities_2d, start_date = calculate_moon_visibility(lon_vals, lat_vals, day, amount, visibility_criterion)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    print_ts(f"Rearranging the visibility values...")
    t1 = time()
    visibilities_3d = np.array(visibilities_2d).reshape(ny, nx, amount)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    print_ts(f"Getting colours for the categories...")
    t1 = time()
    categories, colors_rgba = get_category_colors(visibility_criterion)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    print_ts(f"Categorizing visibilities...")
    t1 = time()
    category_to_index = {cat: i for i, cat in enumerate(categories.keys())}
    visibilities_mapped = map_visibilities(visibilities_3d, category_to_index, ny, nx, amount)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    print_ts(f"Plotting...")
    t1 = time()
    plot_map(lon_vals, lat_vals, visibilities_mapped, states_clip, places_clip, list(categories.keys()), colors_rgba, start_date, amount, path)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

if __name__ == "__main__":
    today = datetime(2024, 8, 4)
    months: int = 1
    map_coord: tuple[int, int, int, int] = (-170, -40, 15, 61)
    resolution: int = 150
    days: int = 3
    criterion: int = 0 # Either 0 (Odeh, 2006), or 1 (Yallop, 1997)

    sys.stdout = Tee(f"mapper_logs/mapper_{datetime.fromtimestamp(time()).strftime("%Y-%m-%d_%H%M%S")}.log")
    start_time: float = time()
    for month in range(months):
        month_start_time: float = time()
        new_day = today + timedelta(days=np.round(29.5 * month))

        # Islamic Date Formatting
        islamic_date = gregorian_to_hijri(new_day.year, new_day.month, new_day.day)
        islamic_year, islamic_month, islamic_day = islamic_date[0], islamic_date[1], islamic_date[2]
        if islamic_day > 6:
            islamic_month += 1
            if islamic_month > 12:
                islamic_month = 1
                islamic_year += 1
        islamic_month_name = get_islamic_month(islamic_month)

        # Create path
        new_path = f"B:/Personal/New Moon Visibilities/Experiment/"
        if not os.path.exists(new_path):
            print_ts(f"Creating {new_path}...")
            os.makedirs(new_path)

        print_ts(f"===Generating map for {islamic_month_name}, {islamic_year}===")
        main(day=new_day, res=resolution, amount=days, visibility_criterion=criterion, path=new_path)
        print_ts(f"===Map for {islamic_month_name}, {islamic_year} Complete===")
        print_ts(f"Time to generate map for {islamic_month_name}, {islamic_year}: {(time() - month_start_time):.2f}s")

    print_ts(f"Total time taken: {(time() - start_time):.2f}s")