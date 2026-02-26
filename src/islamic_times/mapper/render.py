"""Matplotlib/geopandas rendering for visibility maps."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from textwrap import wrap
from typing import Any

import numpy as np

from ._deps import require_mapper_dependencies
from .config import RenderConfig
from .palette import category_entries


def signed_log_transform(x: np.ndarray, epsilon: float) -> np.ndarray:
    return np.sign(x) * np.log1p(np.abs(x) / epsilon)


def inverse_signed_log_transform(y: float, epsilon: float) -> float:
    return float(np.sign(y) * (np.expm1(np.abs(y)) * epsilon))


def _setup_color_mapping(mode: str, visibilities: np.ndarray, criterion: int) -> tuple[Any, Any, float | None, list[str], dict[str, Any]]:
    require_mapper_dependencies()
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    if mode == "raw":
        mask_valid = (~np.isin(visibilities, [-999, -998])) & (~np.isnan(visibilities))
        valid_data = visibilities[mask_valid]
        if valid_data.size == 0:
            raise ValueError("No valid q_values to display in raw mode.")

        epsilon = max(float(np.percentile(np.abs(valid_data), 50)), 0.1)
        transformed_data = signed_log_transform(valid_data, epsilon=epsilon)
        cmap = plt.get_cmap("viridis")
        norm = mcolors.Normalize(vmin=float(np.min(transformed_data)), vmax=float(np.max(transformed_data)))
        return cmap, norm, epsilon, [], {}

    entries = category_entries(criterion)
    labels = [label for label, _, _ in entries]
    rgba = {}
    colors = []
    for label, hex_color, alpha in entries:
        r, g, b, _ = mcolors.to_rgba(hex_color)
        rgba[label] = (r, g, b, alpha)
        colors.append((r, g, b, alpha))

    cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(len(labels) + 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    return cmap, norm, None, labels, rgba


def _plot_features(ax: Any, states_clip: Any, places_clip: Any, annotate_cities: bool) -> None:
    require_mapper_dependencies()
    from matplotlib.patheffects import Normal, Stroke

    states_clip.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=0.65)
    places_clip.plot(ax=ax, color="violet", markersize=7)
    if not annotate_cities:
        return
    for _, row in places_clip.iterrows():
        text = ax.text(row.geometry.x + 0.05, row.geometry.y + 0.05, row["NAME"], fontsize=12, color="white")
        text.set_path_effects([Stroke(linewidth=2, foreground="black"), Normal()])


def _create_legend(fig: Any, gs: Any, labels: list[str], rgba: dict[str, Any]) -> None:
    require_mapper_dependencies()
    from matplotlib.patches import Rectangle

    legend_ax = fig.add_subplot(gs[:, 1])
    legend_ax.axis("off")
    row_height = 1.05
    for idx, label in enumerate(labels):
        wrapped = "\n".join(wrap(label, width=30))
        legend_ax.add_patch(Rectangle((0, idx * row_height), 1, 1, color=rgba[label]))
        legend_ax.text(3.2, idx * row_height + 0.5, wrapped, fontsize=12, va="center", ha="left")
    legend_ax.set_xlim(0, 2)
    legend_ax.set_ylim(0, len(labels) * row_height)


def _create_scale(fig: Any, mesh: Any, norm: Any, epsilon: float) -> None:
    require_mapper_dependencies()
    from matplotlib.patches import Rectangle

    cbar_ax = fig.add_axes([0.88, 0.135, 0.02, 0.8])
    cbar = fig.colorbar(mesh, cax=cbar_ax)
    cbar.set_label("Q Value", fontsize=12)
    ticks = np.linspace(norm.vmin, norm.vmax, 7)
    labels = [f"{inverse_signed_log_transform(float(t), epsilon=epsilon):.1f}" for t in ticks]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)

    legend_ax = fig.add_axes([0.82, 0.04, 0.1, 0.08])
    legend_ax.axis("off")
    for i, (label_text, color) in enumerate(
        [("Moonset before the new moon.", "#141414"), ("Moonset before sunset.", "#393a3c")]
    ):
        y = 1 - i * 0.5
        legend_ax.add_patch(Rectangle((0, y - 0.3), 0.3, 0.3, facecolor=color, edgecolor="black", linewidth=1.2))
        legend_ax.text(0.4, y - 0.15, label_text, fontsize=10, va="center", ha="left", color="black")


def _name_figure(start_date: datetime, islamic_month_name: str, islamic_year: int, criterion: int, mode: str) -> tuple[str, int]:
    name = f"{start_date.strftime('%Y-%m-%d')} {islamic_month_name} {islamic_year}"
    name += "—Yallop" if criterion == 1 else "—Odeh"
    quality = 95 if mode == "raw" else 90
    if mode == "raw":
        name += " Gradient"
    name += ".jpg"
    return name, quality


def render_visibility_map(
    lon_edges: np.ndarray,
    lat_edges: np.ndarray,
    lon_centers: np.ndarray,
    lat_centers: np.ndarray,
    visibilities: np.ndarray,
    states_clip: Any,
    places_clip: Any,
    start_date: datetime,
    islamic_month_name: str,
    islamic_year: int,
    criterion: int,
    render_cfg: RenderConfig,
    out_dir: Path,
) -> Path:
    """Render and write a multi-day map figure to disk."""
    require_mapper_dependencies()
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt

    amount = visibilities.shape[2]
    cmap, norm, epsilon, labels, rgba = _setup_color_mapping(render_cfg.map_mode, visibilities, criterion)

    fig = plt.figure(figsize=(20, 15), dpi=render_cfg.dpi, constrained_layout=False)
    gs = gridspec.GridSpec(amount, 2, width_ratios=[50, 1], height_ratios=[2] * amount)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(amount)]
    mesh = None

    for day_idx, ax in enumerate(axes):
        if render_cfg.map_mode == "raw":
            z_data_raw = visibilities[:, :, day_idx]
            special_mask = np.isin(z_data_raw, [-999, -998])
            z_data = np.where(special_mask, np.nan, z_data_raw)
            assert epsilon is not None
            z_data_transformed = signed_log_transform(z_data, epsilon=epsilon)
            mesh = ax.pcolormesh(lon_edges, lat_edges, z_data_transformed, cmap=cmap, norm=norm, shading="auto", antialiased=False)
        else:
            day_data = visibilities[:, :, day_idx]
            mesh = ax.pcolormesh(lon_edges, lat_edges, day_data, cmap=cmap, norm=norm, shading="auto", antialiased=False)

        _plot_features(ax, states_clip, places_clip, render_cfg.annotate_cities)
        ax.set_xlim(float(np.min(lon_edges)), float(np.max(lon_edges)))
        ax.set_ylim(float(np.min(lat_edges)), float(np.max(lat_edges)))
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"New Moon Visibility on {(start_date + timedelta(days=day_idx)).strftime('%Y-%m-%d')} at Local Best Time")
        ax.set_xticks(lon_centers[:: max(1, len(lon_centers) // 10)])
        ax.set_yticks(lat_centers[:: max(1, len(lat_centers) // 10)])

    if render_cfg.map_mode == "category":
        _create_legend(fig, gs, labels, rgba)
    else:
        assert mesh is not None and epsilon is not None
        _create_scale(fig, mesh, norm, epsilon)

    criterion_string = "Odeh, 2006" if criterion == 0 else "Yallop, 1997"
    plt.subplots_adjust(hspace=0.2, left=0.05, right=0.85, top=0.95, bottom=0.05)
    plt.figtext(0.15, 0.01, f"The New Moon (i.e. conjunction) occurs at {start_date.strftime('%Y-%m-%d %X')} UTC", ha="center", fontsize=12)
    plt.figtext(0.945, 0.03, f"Criterion: {criterion_string}", ha="center", fontsize=12)
    plt.figtext(0.84, 0.01, "CC BY-SA | Hassan Tahan | Created with the islamic_times Python library", ha="center", fontsize=12)
    plt.figtext(
        0.5,
        0.98,
        f"{amount}-Day New Moon Crescent Visibility Map for {islamic_month_name}, {islamic_year} A.H.",
        ha="center",
        fontsize=16,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    filename, default_quality = _name_figure(start_date, islamic_month_name, islamic_year, criterion, render_cfg.map_mode)
    quality = render_cfg.jpeg_quality or default_quality
    output_path = out_dir / filename
    plt.savefig(output_path, format=render_cfg.image_format, pil_kwargs={"optimize": True, "progressive": True, "quality": quality})
    plt.close(fig)
    return output_path

