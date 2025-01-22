from datetime import timedelta, datetime
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import matplotlib.patheffects as path_effects
from shapely.geometry import Polygon
from islamic_times import islamic_times as it

cities = [
        'Washington,  D.C.', 'Los Angeles', 'London',
        'Vienna', 'Moscow', 'Tokyo', 'Beijing', 'Mumbai',
        'Cairo', 'Lagos', 'Sydney', 'São Paulo', 'Mexico City', 
        'Toronto', 'Istanbul', 'Tehran', 'Islamabad', 'Perth',
        'Bangkok', 'Hong Kong', 'Singapore', 'Makkah', 'Miami',
        'Lima', 'Bogota', 'Cape Town', 'Madrid', 'Vancouver',
        'Nairobi', 'Addis Ababa', 'Dakar', 'Santiago'
    ]

def compute_q_values_for_lat(i_lat: float, lon: float, day: datetime, amount: int, type: int):
        q_values_temp = []
        for i_lon in lon:
            local = it.ITLocation(latitude=i_lat, longitude=i_lon, elevation=0, today=day, find_local_tz=False)
            vis = local.visibilities(days=amount, type=type)
            for label, value in vis.items():
                q_values_temp.append(value[1])
        return q_values_temp

def calculate_moon_visibility(lon: float, lat: float, day: datetime, amount: int = 1, type: int = 0):
    moon_phases = it.ITLocation(latitude=0, longitude=0, today=day, find_local_tz=False).moonphases()
    for item in moon_phases:
            if item['phase'] == "New Moon":
                new_moon_date = item['datetime']  
    
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

def main(day: datetime, res: int, amount: int = 1, type: int = 0):
    # Load shapefiles for states/provinces and populated places
    states_provinces_path = "map_shp_files/combined_polygons.shp"
    places_path = "map_shp_files/combined_points.shp"

    states_gdf = gpd.read_file(states_provinces_path)

    # Filter for listed cities
    places_gdf = gpd.read_file(places_path)
    places_gdf = places_gdf[places_gdf['NAME'].isin(cities)]
    places_gdf = (
        places_gdf.loc[places_gdf.groupby('NAME')['POP_MAX'].idxmax()]
    )

    # Longitudes
    minx, maxx = -179, 180 
    # Latitudes
    miny, maxy = -61, 61

    # Clip map to lat & long
    bbox_polygon = Polygon([
        (minx, miny),
        (maxx, miny),
        (maxx, maxy),
        (minx, maxy)
    ])
    bbox_gdf = gpd.GeoDataFrame([1], geometry=[bbox_polygon], crs=states_gdf.crs)
    states_clip = gpd.overlay(states_gdf, bbox_gdf, how='intersection')
    places_clip = places_gdf.cx[minx:maxx, miny:maxy]

    # Create grid
    resolution = res
    lon_vals = np.linspace(minx, maxx, resolution)
    lat_vals = np.linspace(miny, maxy, resolution)
    nx = len(lon_vals)
    ny = len(lat_vals)

    # Compute the visibilities
    visibilities_2d, start_date = calculate_moon_visibility(
        lon_vals, lat_vals, day, amount=amount, type=type
    )
    # Should have shape (ny, nx*amount).
    visibilities_2d = np.array(visibilities_2d)

    # Reshape to 3D: (ny, nx, amount)
    # Each row has (nx * amount) columns, split into 'nx' chunks,
    # and in each chunk, we have 'amount' classification strings.
    if visibilities_2d.shape == (ny, nx * amount):
        # reshape each row from (nx * amount) to (nx, amount)
        # then stack those along the lat dimension => (ny, nx, amount)
        visibilities_3d = visibilities_2d.reshape(ny, nx, amount)
    else:
        raise ValueError(f"Unexpected shape {visibilities_2d.shape} - "
                        f"expected ({ny}, {nx*amount}).")
    
    # Category colour mapping
    if type == 0:
        category_colors = {
            "Moonset before the new moon.": "#1A1A1A",
            "Moonset before sunset.": "#303030",
            "D: Crescent is not visible even by optical aid.": "#807f80",
            "C: Crescent is visible by optical aid only.": "#B89D18",
            "B: Crescent is visible by optical aid, and it could be seen by naked eyes.": "#74b818",
            "A: Crescent is visible by naked eyes.": "#1BB818"
        }
        transparent_category = "D: Crescent is not visible even by optical aid."
    else:
        category_colors = {
            "Moonset before the new moon.": "#1A1A1A",
            "Moonset before sunset.": "#303030",
            "F: Not visible; below the Danjon limit.": "#807f80",
            "E: Not visible with a [conventional] telescope.": "#B81818",
            "D: Will need optical aid to find crescent.": "#e3d61b",
            "C: May need optical aid to find crescent.": "#89d518",
            "B: Visible under perfect conditions.": "#54b818",
            "A: Easily visible.": "#1bdf18",
        }
        transparent_category = "F: Not visible; below the Danjon limit."

    category_colors_rgba = {
        cat: mcolors.to_rgba(color) if cat != transparent_category else (0, 0, 0, 0.1)
        for cat, color in category_colors.items()
    }
    unique_categories = list(category_colors.keys())
    category_to_index = {cat: i for i, cat in enumerate(unique_categories)}

    # Convert from string categories to integer indices
    visibilities_mapped = np.zeros((ny, nx, amount), dtype=int)

    for i_lat in range(ny):
        for i_lon in range(nx):
            for i_day in range(amount):
                cat_str = visibilities_3d[i_lat, i_lon, i_day]
                if cat_str not in category_to_index:
                    raise ValueError(f"Unknown category '{cat_str}'")
                visibilities_mapped[i_lat, i_lon, i_day] = category_to_index[cat_str]

    cmap = mcolors.ListedColormap([category_colors_rgba[cat] for cat in unique_categories])
    bounds = np.arange(len(unique_categories) + 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Create figure and grid layout
    fig = plt.figure(figsize=(20, 15), dpi=80)
    gs = gridspec.GridSpec(amount, 2, width_ratios=[50, 1], height_ratios=[2] * amount)  # Adjust grid layout

    # Create subplots for maps
    axes = [fig.add_subplot(gs[i, 0]) for i in range(amount)]

    for i_day, ax in enumerate(axes):
        # Slice out shape (ny, nx) => lat-by-lon
        day_slice = visibilities_mapped[..., i_day]

        # Plot each day with pcolormesh
        mesh = ax.pcolormesh(
            lon_vals,
            lat_vals,
            day_slice,
            cmap=cmap,
            norm=norm,
            shading="auto"
        )

        # Plot states/provinces
        states_clip.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=0.5)

        # Plot places
        places_clip.plot(ax=ax, color='violet', markersize=7)

        # Label cities with black outline and white filling
        for idx, row in places_clip.iterrows():
            text = ax.text(
                row.geometry.x + 0.05,  # Offset slightly to the right
                row.geometry.y + 0.05,  # Offset slightly upward
                row['NAME'],            # City name column
                fontsize=12,
                color='white',          # Text filling color
                ha='left',
                va='bottom'
            )

            # Add black outline to the text
            text.set_path_effects([
                path_effects.Stroke(linewidth=2, foreground='black'),
                path_effects.Normal()
            ])

        ax.set_xlim(minx, maxx)
        ax.set_ylim(miny, maxy)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

        # Add title for each subplot
        day_date = start_date + timedelta(days=i_day)
        ax.set_title(f"New Moon Visibility on {day_date.strftime('%d-%m-%Y')} at Local Best Time.")

    # Add a legend with color boxes to the right
    legend_ax = fig.add_subplot(gs[:, 1])
    legend_ax.axis('off')  # Turn off the axis

    # Create legend items
    for idx, category in enumerate(unique_categories):
        def wrap_text(text, width=10):
            import textwrap
            return "\n".join(textwrap.wrap(text, width))
        
        wrapped_category = wrap_text(category, width=15)
        legend_ax.add_patch(
            Rectangle((0, idx), 1, 1, color=category_colors_rgba[category])
        )
        legend_ax.text(
            1.2, idx + 0.5, wrapped_category, fontsize=12, va='center', ha='left'
        )

    legend_ax.set_xlim(0, 2)
    legend_ax.set_ylim(0, len(unique_categories))
    legend_ax.set_aspect('auto')

    plt.subplots_adjust(hspace=0.2)
    plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.05)

    # Add caption below the maps
    new_moon_date_utc = start_date.strftime("%H:%M:%S %Y-%m-%d")
    fig.text(0.5, 0.01, f"Date of New Moon at UTC: {new_moon_date_utc}",
            ha="center", fontsize=16)

    plt.show()


if __name__ == "__main__":
    main(datetime(2025, 6, 25), res=200, amount=3, type=1)