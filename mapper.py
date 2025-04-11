import os
import sys
import numpy as np
from time import time
import geopandas as gpd
from typing import List, Tuple
from multiprocessing import Pool, cpu_count
from datetime import timedelta, datetime
import islamic_times.astro_core as fast_astro
from islamic_times.time_equations import get_islamic_month, gregorian_to_hijri

# Plotting libraries
from textwrap import wrap
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
from matplotlib.patheffects import Stroke, Normal
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

AVERAGE_LUNAR_MONTH_DAYS: int = 29.53059

CITIES_WORLD: list[str] = [
        # PACIFIC
        'Honolulu', 
        
        # NORTH AMERICA
        'Vancouver', 'Los Angeles', 'Mexico City',
        'Toronto', 'Miami', 'Washington,  D.C.',

        # SOUTH AMERICA
        'Lima', 'Bogota', 'Santiago', 'São Paulo',

        # WEST AFRICA
        'Dakar', 'Lagos',

        # EUROPE
        'Madrid', 'London', 'Vienna', 'Moscow',

        # SOUTH AFRICA
        'Cape Town', 

        # MIDDLE EAST
        'Istanbul', 'Cairo', 'Makkah', 'Tehran',

        # EAST AFRICA
        'Nairobi', 'Addis Ababa', 

        # SOUTH ASIA
        'Islamabad', 'Mumbai',

        # SOUTH EAST ASIA
        'Bangkok',  'Singapore',

        # EAST ASIA
        'Hong Kong', 'Beijing', 'Tokyo', 
        
        # AUSTRALIA
        'Sydney', 'Perth' 
    ]

CITIES_IRAN: list[str] = [
    'Tehran', 'Mashhad', 'Kerman', 'Shiraz', 'Zanjan', 
    'Ardabil', 'Isfahan', 'Gorgan', 'Tabriz', 'Semnan', 
    'Yazd', 'Rasht', 'Arak', 'Boshruyeh', 'Mehran', 'Dargaz', 
    'Chabahar', 'Zahedan', 'Birjand', 'Sanandaj', 'Ahvaz', 
    'Saravan', 'Hamadan', 'Khorramabad', 'Qomsheh', 'Ilam', 
    'Sari', 'Qazvin', 'Bandar-e-Abbas', 'Bandar-e Bushehr', 
    'Sirjan', 'Kashmar', 'Bojnurd', 'Qom', 'Urmia', 'Khvoy',
    'Yasuj'
    ]

CITIES_MIDDLE_EAST: list[str] = [
    'Istanbul', 'Khartoum', 'Cairo', 'Luxor', 'Ankara', 
    'Beirut', 'Aleppo', 'Medina', 'Makkah', 'Djibouti',
    'Sanaa', 'Irbil', 'Baghdad', 'Riyadh', 'Kuwait City',
    'Baku', 'Tehran', 'Doha', 'Dubai', 'Kerman', 'Muscat', 
    'Mashhad', 'Karachi', 'Kabul'
    ]

CITIES_NORTH_AMERICA: list[str] = [
    # Pacific
    'Honolulu',
    
    # CANADA
    'Vancouver', 'Edmonton', 'Calgary', 'Winnipeg', 'Thunder Bay', 
    'Toronto', 'Montréal', 'Halifax', 'St. John\'s',

    # UNITED STATES
    'Portland', 'San Francisco', 'Los Angeles', 'Billings', 
    'Albuquerque', 'Denver', 'Kansas City', 'Dallas', 'Houston', 
    'Minneapolis', 'Chicago', 'Orlando', 'Atlanta', 'Miami', 
    'Washington,  D.C.', 'Boston',

    # MEXICO
    'Hermosillo', 'Monterrey', 'Mexico City', 'Mérida',

    # CARIBBEAN
    'Havana', 'Kingston',
]

CITIES_EUROPE: list[str] = [
    'Lisbon', 'Dublin', 'Madrid', 'Edinburgh',
    'London', 'Barcelona', 'Paris', 'Amsterdam', 'Zürich', 
    'Oslo', 'Rome', 'København', 'Venice', 'Berlin', 
    'Vienna', 'Stockholm', 'Sarajevo', 'Warsaw', 'Athens',
    'Riga', 'Bucharest', 'Minsk', 'Istanbul', 'Kyiv',
    'Ankara', 'Moscow', 'Rostov', 'Tbilisi'
]

REGION_COORDINATES: dict[str, tuple[int, int, int, int]] = {
    'WORLD'         :   (-179, 180, -61, 61),
    'WORLD_FULL'    :   (-179, 180, -89, 90),
    'NORTH_AMERICA' :   (-170, -40, 15, 61),
    'EUROPE'        :   (-15, 50, 34, 61),
    'MIDDLE_EAST'   :   (25, 75, 10, 45),
    'IRAN'          :   (43.5, 63.5, 24.5, 40)
}

REGION_CITIES: dict[str, list[str]] = {
    'WORLD'         : CITIES_WORLD,
    'WORLD_FULL'    : CITIES_WORLD, # NOT YET SUPPORTED
    'NORTH_AMERICA' : CITIES_NORTH_AMERICA,
    'EUROPE'        : CITIES_EUROPE,
    'MIDDLE_EAST'   : CITIES_MIDDLE_EAST,
    'IRAN'          : CITIES_IRAN 
}

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

def print_ts(message: str):
    print(f"[{datetime.fromtimestamp(time()).strftime('%X %d-%m-%Y')}] {message}")

def split_lat_chunks(lat_vals, n_chunks):
    return np.array_split(lat_vals, n_chunks)

def batch_worker(lat_chunk, lon_vals, dt, days, criterion, utc_offset, elev, temp, press, mode_byte):
    lats, lons = np.meshgrid(lat_chunk, lon_vals, indexing="ij")
    lats_flat = lats.ravel()
    lons_flat = lons.ravel()
    result_flat = fast_astro.compute_visibilities_batch(lats_flat, lons_flat, dt, days, criterion,
                                                utc_offset, elev, temp, press, mode_byte)
    ny, nx = lats.shape
    return result_flat.reshape(ny, nx, days)

def compute_visibility_map_parallel(lon_vals, lat_vals, new_moon_date, days, criterion,
                                    utc_offset=0.0, elev=0.0, temp=20.0, press=101.325,
                                    mode="category", num_workers=cpu_count()):
    mode_byte: str = 'r' if mode == "raw" else 'c'
    lat_chunks = split_lat_chunks(lat_vals, num_workers)

    print_ts(f"Conjunction Date: {new_moon_date.strftime("%Y-%m-%d %X")}")

    args = [(chunk, lon_vals, new_moon_date, days, criterion, utc_offset, elev, temp, press, mode_byte)
            for chunk in lat_chunks]

    with Pool(num_workers) as pool:
        results = pool.starmap(batch_worker, args)

    # Stack along latitude axis to reconstruct full (ny, nx, days) array
    visibilities_3d = np.vstack(results)
    return visibilities_3d

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

def map_visibilities(visibilities_3d, category_to_index, ny, nx, amount, mode="category"):
    if mode == "raw":
        return visibilities_3d  # already numeric
    else:
        visibilities_mapped = np.zeros((ny, nx, amount), dtype=int)
        for i_lat in range(ny):
            for i_lon in range(nx):
                for i_day in range(amount):
                    visibilities_mapped[i_lat, i_lon, i_day] = category_to_index[visibilities_3d[i_lat, i_lon, i_day]]
        return visibilities_mapped

def signed_log_transform(x, epsilon: float):
    """Apply a signed pseudo-log transform to handle both negative and positive values."""
    return np.sign(x) * np.log1p(np.abs(x) / epsilon)

def inverse_signed_log_transform(y, epsilon: float):
    """Reverse the signed_log_transform to recover original q_value from transformed."""
    return np.sign(y) * (np.expm1(np.abs(y)) * epsilon)

def setup_color_mapping(mode, visibilities_mapped, unique_categories, category_colors_rgba):
    if mode == "raw":
        # Filter all valid values once across all days
        mask_valid = (~np.isin(visibilities_mapped, [-999, -998])) & (~np.isnan(visibilities_mapped))
        valid_data = visibilities_mapped[mask_valid]

        if valid_data.size == 0:
            raise ValueError("No valid q_values to display in raw mode.")

        # Use median of abs(q) values for epsilon
        epsilon = np.percentile(np.abs(valid_data), 50)
        epsilon = max(epsilon, 0.1)

        # Apply signed log transform to all valid values
        transformed_data = signed_log_transform(valid_data, epsilon=epsilon)

        # Global min and max for color normalization
        zmin_transformed = np.min(transformed_data)
        zmax_transformed = np.max(transformed_data)

        cmap = plt.get_cmap("viridis")
        norm = mcolors.Normalize(vmin=zmin_transformed, vmax=zmax_transformed)

        return cmap, norm, epsilon
    else:
        cmap = mcolors.ListedColormap([category_colors_rgba[cat] for cat in unique_categories])
        bounds = np.arange(len(unique_categories) + 1)
        norm = mcolors.BoundaryNorm(bounds, cmap.N)
        return cmap, norm, None

def plot_raw_map(ax, lon_vals, lat_vals, z_data_raw, cmap, epsilon, norm):
    # Create a mask for the special-case values only once.
    special_mask = np.isin(z_data_raw, [-999, -998])
    z_data = np.where(special_mask, np.nan, z_data_raw)

    # Transform valid q_values using the provided epsilon.
    z_data_transformed = signed_log_transform(z_data, epsilon=epsilon)

    # Use shared/global norm from setup_color_mapping()
    mesh = ax.pcolormesh(lon_vals, lat_vals, z_data_transformed, cmap=cmap, norm=norm, shading="auto")

    # Fill special-case areas once per value.
    for value, color in [(-999, '#141414'), (-998, '#393a3c')]:
        mask = (z_data_raw == value).astype(float)
        if np.any(mask):
            ax.contourf(lon_vals, lat_vals, mask, levels=[0.5, 1.5], colors=[color], alpha=1.0)

    # Draw contours using same norm
    valid = z_data_transformed[~np.isnan(z_data_transformed)]
    if valid.size > 0:
        contour_levels = np.linspace(norm.vmin, norm.vmax, 10)
        cs = ax.contour(lon_vals, lat_vals, z_data_transformed, levels=contour_levels,
                        colors='white', linewidths=1.2)
        fmt = {lvl: f"{inverse_signed_log_transform(lvl, epsilon=epsilon):.1f}" for lvl in cs.levels}
        ax.clabel(cs, cs.levels, fmt=fmt, inline=True, fontsize=10)

    return mesh

def plot_features(ax, states_clip, places_clip):
    states_clip.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=0.65)
    places_clip.plot(ax=ax, color='violet', markersize=7)
    # For many points, consider vectorizing or batching text annotations if performance remains an issue.
    for _, row in places_clip.iterrows():
        text = ax.text(row.geometry.x + 0.05, row.geometry.y + 0.05, row['NAME'],
                       fontsize=12, color='white', ha='left', va='bottom')
        text.set_path_effects([Stroke(linewidth=2, foreground='black'), Normal()])

def create_legend(fig, gs, unique_categories, category_colors_rgba):
    legend_ax = fig.add_subplot(gs[:, 1])
    legend_ax.axis("off")

    def wrap_text(text, width=20):
        return "\n".join(wrap(text, width))

    row_height = 1.05
    for idx, category in enumerate(unique_categories):
        wrapped = wrap_text(category, width=30)
        legend_ax.add_patch(Rectangle((0, idx * row_height), 1, 1, color=category_colors_rgba[category]))
        legend_ax.text(3.2, idx * row_height + 0.5, wrapped, fontsize=12, va='center', ha='left')

    legend_ax.set_xlim(0, 2)
    legend_ax.set_ylim(0, len(unique_categories) * row_height)
    legend_ax.set_aspect("auto")

def create_scale(fig, mesh, norm):
    cbar_ax = fig.add_axes([0.88, 0.135, 0.02, 0.8])
    cbar = fig.colorbar(mesh, cax=cbar_ax)
    cbar.set_label("Q Value", fontsize=12)

    tick_values_transformed = np.linspace(norm.vmin, norm.vmax, 7)
    tick_labels = [f"{inverse_signed_log_transform(t, epsilon=2.0):.1f}" for t in tick_values_transformed]
    cbar.set_ticks(tick_values_transformed)
    cbar.set_ticklabels(tick_labels)

    # Special cases
    legend_ax: Axes = fig.add_axes([0.82, 0.04, 0.1, 0.08])  # [left, bottom, width, height]
    legend_ax.axis("off")
    special_cases = [
        ("Moonset before the new moon.", "#141414"),
        ("Moonset before sunset.", "#393a3c"),
    ]
    for i, (label_text, color) in enumerate(special_cases):
        y = 1 - i * 0.5
        legend_ax.add_patch(Rectangle((0, y - 0.3), 0.3, 0.3,
                                      facecolor=color,
                                      edgecolor="black",
                                      linewidth=1.2))
        legend_ax.text(0.4, y - 0.15, label_text, fontsize=10, va='center', ha='left', color='black')

def annotate_plot(fig, start_date, criterion, days, islamic_month_name, islamic_year):
    criterion_string = "Odeh, 2006" if criterion == 0 else "Yallop, 1997"
    plt.subplots_adjust(hspace=0.2, left=0.05, right=0.85, top=0.95, bottom=0.05)
    plt.figtext(0.15, 0.01,
                f"The New Moon (i.e. conjunction) occurs at {start_date.strftime('%Y-%m-%d %X')} UTC",
                ha="center", fontsize=12)
    plt.figtext(0.945, 0.03, f"Criterion: {criterion_string}", ha="center", fontsize=12)
    plt.figtext(0.84, 0.01,
                "CC BY-SA | Hassan Tahan | Created with the islamic_times Python library",
                ha="center", fontsize=12)
    plt.figtext(0.5, 0.98,
                f"{days}-Day New Moon Crescent Visibility Map for {islamic_month_name}, {islamic_year} A.H.",
                ha="center", fontsize=16)

def name_fig(start_date, islamic_month_name, islamic_year, mode):
    name = f"{start_date.strftime('%Y-%m-%d')} {islamic_month_name} {islamic_year}"
    name += "—Yallop" if criterion == 1 else "Odeh"
    qual = 95 if mode == "raw" else 90
    if mode == "raw":
        name += " Gradient"
    name += ".jpg"
    return name, qual

def plot_map(lon_vals, lat_vals, visibilities_mapped, states_clip, places_clip,
             unique_categories, category_colors_rgba, start_date, amount, out_dir, 
             islamic_month_name, islamic_year, mode="category"):
    # Set up the color mapping and obtain epsilon if in raw mode.
    print_ts("Plotting: Setting up colour map...")
    cmap, norm, epsilon = setup_color_mapping(mode, visibilities_mapped, unique_categories, category_colors_rgba)

    print_ts("Plotting: Adding subplots...")
    width_x, width_y = 20, 15
    dpi = 300
    fig = plt.figure(figsize=(width_x, width_y), dpi=dpi, constrained_layout=False)
    gs = gridspec.GridSpec(amount, 2, width_ratios=[50, 1], height_ratios=[2] * amount)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(amount)]
    mesh = None

    # Plot each day's visibility.
    for i_day, ax in enumerate(axes):
        print_ts(f"Plotting: Plotting Day {i_day + 1} ...")
        if mode == "raw":
            z_data_raw = visibilities_mapped[..., i_day]
            print_ts(f"Plotting: Raw map plotting for Day {i_day + 1} ...")
            mesh = plot_raw_map(ax, lon_vals, lat_vals, z_data_raw, cmap, epsilon, norm)
        else:
            data = visibilities_mapped[..., i_day]
            mesh = ax.pcolormesh(lon_vals, lat_vals, data, cmap=cmap, norm=norm, shading="auto")

        print_ts(f"Plotting: Adding features for {i_day + 1} ...")
        plot_features(ax, states_clip, places_clip)
        ax.set_xlim(min(lon_vals), max(lon_vals))
        ax.set_ylim(min(lat_vals), max(lat_vals))
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"New Moon Visibility on {(start_date + timedelta(days=i_day)).strftime('%Y-%m-%d')} at Local Best Time")

    if mode == "category":
        print_ts("Plotting: Adding legent...")
        create_legend(fig, gs, unique_categories, category_colors_rgba)
    else:
        print_ts("Plotting: Adding scale...")
        create_scale(fig, mesh, norm)

    print_ts("Plotting: Annotating plot...")
    annotate_plot(fig, start_date, criterion, days_to_generate, islamic_month_name, islamic_year)
    name, qual = name_fig(start_date, islamic_month_name, islamic_year, mode)

    print_ts("Plotting: Saving...")
    plt.savefig(os.path.join(out_dir, name), format='jpg',
                pil_kwargs={'optimize': True, 'progressive': True, 'quality': qual})
    plt.close()

def plotting_loop(new_moon_date: datetime, map_params: Tuple, mode: str = "category", region: str = 'WORLD', 
                  amount: int = 1, visibility_criterion: int = 0):
    # Start timing for the month
    month_start_time: float = time()
    
    # Islamic Date Formatting
    islamic_date = gregorian_to_hijri(new_moon_date.year, new_moon_date.month, new_moon_date.day)
    islamic_year, islamic_month, islamic_day = islamic_date[0], islamic_date[1], islamic_date[2]
    if islamic_day > 6:
        islamic_month += 1
        if islamic_month > 12:
            islamic_month = 1
            islamic_year += 1
    islamic_month_name = get_islamic_month(islamic_month)

    # Unpack map parameters
    states_clip, places_clip, lon_vals, lat_vals, nx, ny = map_params

    # Create path
    path = f"{master_path}{map_region.replace('_', ' ').title()}/{islamic_year}/"
    if not os.path.exists(path):
        print_ts(f"Creating {path}...")
        os.makedirs(path)

    # Start
    print_ts(f"===Generating map for {islamic_month_name}, {islamic_year}===")

    # Calculate
    print_ts(f"Calculating new moon crescent visibilities...")
    t1 = time()
    visibilities_3d = compute_visibility_map_parallel(lon_vals, lat_vals, new_moon_date, amount, visibility_criterion, mode=mode)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    # Categorization
    print_ts(f"Getting colours for the categories...")
    t1 = time()
    categories, colors_rgba = get_category_colors(visibility_criterion)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    # Categorize
    print_ts(f"Colouring visibilities...")
    t1 = time()
    category_to_index = {cat: i for i, cat in enumerate(categories.keys())}
    visibilities_mapped = map_visibilities(visibilities_3d, category_to_index, ny, nx, amount, mode=mode)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    # Plotting
    print_ts(f"Plotting...")
    t1 = time()
    plot_map(lon_vals, lat_vals, visibilities_mapped, states_clip, places_clip, 
             list(categories.keys()) if mode == "category" else [], 
             colors_rgba if mode == "category" else {}, 
             new_moon_date, amount, path, 
             islamic_month_name, islamic_year, mode)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    # Finished
    print_ts(f"===Map for {islamic_month_name}, {islamic_year} Complete===")
    print_ts(f"Time to generate map for {islamic_month_name}, {islamic_year}: {(time() - month_start_time):.2f}s")

def main():
    sys.stdout = Tee(f'mapper_logs/mapper_{datetime.fromtimestamp(time()).strftime("%Y-%m-%d_%H%M%S")}.log')
    start_time: float = time()

    # Select region 
    cities: List[str] = REGION_CITIES[map_region]
    coords: Tuple[int, int, int, int] = REGION_COORDINATES[map_region]

    # Generate the map from the shp files
    print_ts(f"Loading shape files...")
    t1 = time()
    states_gdf, places_gdf = load_shapefiles("map_shp_files/combined_polygons.shp", "map_shp_files/combined_points.shp", cities)
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    # Clip
    print_ts(f"Clipping map to coordinates...")
    t1 = time()
    states_clip, places_clip = clip_map(states_gdf, places_gdf, minx=coords[0], maxx=coords[1], miny=coords[2], maxy=coords[3])
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    print_ts(f"Creating map grid...")
    t1 = time()
    lon_vals, lat_vals, nx, ny = create_grid(resolution, minx=coords[0], maxx=coords[1], miny=coords[2], maxy=coords[3])
    print_ts(f"Time taken: {(time() - t1):.2f}s")

    for month in range(total_months):
        assert map_mode in ("raw", "category"), f"Invalid mode: {map_mode}"
        new_moon_date: datetime = fast_astro.next_phases_of_moon_utc(today + timedelta(days=month * AVERAGE_LUNAR_MONTH_DAYS))[0]

        plotting_loop(new_moon_date, map_params=(states_clip, places_clip, lon_vals, lat_vals, nx, ny), region=map_region, amount=days_to_generate, 
                      visibility_criterion=criterion, mode=map_mode)

    print_ts(f"~~~ --- === Total time taken: {(time() - start_time):.2f}s === --- ~~~")

if __name__ == "__main__":
    # Vars
    today = datetime(1996, 9, 12) # .now() - timedelta(days=365.25*30)
    master_path: str = "B:/Personal/New Moon Visibilities/Experiment/C-Rewrite/Full-Test/"
    total_months: int = 1
    map_region: str = "WORLD" # 'NORTH_AMERICA' 'EUROPE' 'MIDDLE_EAST' 'IRAN' 'WORLD' 'WORLD_FULL'
    map_mode: str = "category" #"raw" # "category"
    resolution: int = 300
    days_to_generate: int = 3
    criterion: int = 1 # Either 0 (Odeh, 2006), or 1 (Yallop, 1997)

    map_region = map_region.upper()
    main()