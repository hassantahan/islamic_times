"""CLI for mapper generation and benchmarking."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone

from .config import ComputeConfig, MapperConfig, RenderConfig
from .pipeline import (
    build_visibility_cache,
    generate_maps,
    load_visibility_cache_tabular,
    render_maps_from_cache,
    save_visibility_cache_tabular,
)
from .. import prayer_times as pt


_CRITERION_INFO: dict[int, tuple[str, str]] = {
    0: (
        "Odeh (2006)",
        "Topocentric geometry with Odeh visibility thresholds.",
    ),
    1: (
        "Yallop (1997)",
        "Geocentric q-values with HMNAO TN 69 classification bands.",
    ),
    2: (
        "Shaukat (n.d.)",
        "Yallop q-values with Shaukat-specific classification thresholds.",
    ),
}


class MapperArgumentParser(argparse.ArgumentParser):
    """ArgumentParser with mapper-specific actionable validation hints."""

    def error(self, message: str) -> None:
        lower = message.lower()
        if "argument --criterion" in lower:
            message += (
                "\nHint: supported values are 0 (Odeh), 1 (Yallop), and 2 (Shaukat). "
                "Run --list-criteria for details."
            )
        elif "argument --map_mode" in lower:
            message += "\nHint: --map_mode must be one of: raw, category."
        super().error(message)


def _format_method_listing() -> str:
    rows: list[str] = ["Available prayer methods (primary key -> display name):"]
    for method in pt.prayer_method_catalog():
        primary_key = method.keys[0] if method.keys else method.name.upper()
        rows.append(f"{primary_key:<12} {method.name}")
    return "\n".join(rows)


def _format_criteria_listing() -> str:
    rows: list[str] = ["Available visibility criteria:"]
    for code in sorted(_CRITERION_INFO):
        name, _ = _CRITERION_INFO[code]
        rows.append(f"{code}: {name}")
    return "\n".join(rows)


def _format_criterion_explanation(criterion: int) -> str:
    name, detail = _CRITERION_INFO[criterion]
    return f"Criterion {criterion}: {name}\n{detail}"


def _build_parser() -> argparse.ArgumentParser:
    parser = MapperArgumentParser(
        description="Generate new-moon visibility maps", argument_default=argparse.SUPPRESS,
    )
    subparsers = parser.add_subparsers(dest="command", required=False)

    def add_common_args(target: argparse.ArgumentParser, include_perf_report: bool = True) -> None:
        target.add_argument("--date", help="ISO reference date for choosing the first lunation")
        target.add_argument("--conjunction", help="Exact first conjunction as an ISO datetime with UTC offset")
        target.add_argument("--master_path", "--master-path", help="Directory where maps are written")
        target.add_argument("--total_months", "--total-months", type=int)
        target.add_argument("--map_region", "--map-region")
        target.add_argument("--map_mode", "--map-mode", choices=("raw", "category"))
        target.add_argument("--resolution", type=int)
        target.add_argument("--days_to_generate", "--days", type=int)
        target.add_argument("--criterion", type=int, choices=(0, 1, 2))
        target.add_argument(
            "--category-compute", choices=("dense", "adaptive"),
            help="Category grid solver. Dense evaluates every cell, adaptive samples blocks.",
        )
        target.add_argument("--save_logs", "--save-logs", action="store_true")
        target.add_argument("--max_workers", "--max-workers", type=int)
        if include_perf_report:
            target.add_argument("--perf_report", "--perf-report", help="JSON output path for timing report")

    def add_info_args(target: argparse.ArgumentParser) -> None:
        target.add_argument(
            "--list-methods",
            action="store_true",
            help="List available prayer methods and exit.",
        )
        target.add_argument(
            "--list-criteria",
            action="store_true",
            help="List supported visibility criteria and exit.",
        )
        target.add_argument(
            "--explain-criterion",
            type=int,
            choices=(0, 1, 2),
            help="Print details for one visibility criterion and exit.",
        )

    common = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
    cache_common = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
    add_common_args(common)
    add_common_args(cache_common, include_perf_report=False)
    add_common_args(parser)
    add_info_args(common)
    add_info_args(parser)

    subparsers.add_parser("generate", parents=[common], help="Generate mapper outputs")
    subparsers.add_parser("benchmark", parents=[common], help="Generate outputs and emit benchmark JSON")
    cache_parser = subparsers.add_parser("cache", parents=[cache_common], help="Compute and save visibility cache")
    cache_parser.add_argument(
        "--cache_path", "--cache-path",
        type=str,
        required=True,
        help="Output cache directory for the tabular cache.",
    )
    export_csv_parser = subparsers.add_parser("export-csv", parents=[cache_common], help="Compute and save CSV visibility data")
    export_csv_parser.add_argument(
        "--csv_path", "--csv-path",
        type=str,
        required=True,
        help="Output directory for the CSV export.",
    )

    render_cache_parser = subparsers.add_parser(
        "render-cache", argument_default=argparse.SUPPRESS,
        help="Render maps from precomputed cache without recomputing visibilities",
    )
    render_cache_parser.add_argument(
        "--cache_path", "--cache-path",
        type=str,
        required=True,
        help="Path to cache directory or manifest.json",
    )
    render_cache_parser.add_argument("--master_path", "--master-path", help="Override the output directory")
    render_cache_parser.add_argument("--save_logs", "--save-logs", action="store_true")
    render_cache_parser.add_argument("--perf_report", "--perf-report", help="JSON output path for timing report")
    # Defaults belong to the root parser so a subcommand does not overwrite earlier flags.
    parser.set_defaults(
        date=None, conjunction=None, master_path=None, total_months=1, map_region="WORLD",
        map_mode="category", resolution=300, days_to_generate=3, criterion=1, category_compute=None,
        save_logs=False, max_workers=None, perf_report=None, list_methods=False, list_criteria=False,
        explain_criterion=None, command=None,
    )
    return parser


def _parse_iso_datetime(value: str) -> datetime:
    """Parse an ISO reference date or timestamp, including a UTC Z suffix on Python 3.10."""
    if value.endswith("Z"):
        value = value[:-1] + "+00:00"
    return datetime.fromisoformat(value)


def _config_from_args(args: argparse.Namespace) -> MapperConfig:
    if args.criterion not in _CRITERION_INFO:
        raise ValueError(
            f"Invalid criterion '{args.criterion}'. Supported values are 0 (Odeh), 1 (Yallop), and 2 (Shaukat)."
        )
    date = _parse_iso_datetime(args.date) if args.date else datetime.now(timezone.utc)
    if date.utcoffset() is not None:
        date = date.astimezone(timezone.utc)
    conjunction = _parse_iso_datetime(args.conjunction) if args.conjunction else None
    default_compute = ComputeConfig()
    return MapperConfig(
        date=date,
        total_months=args.total_months,
        map_region=args.map_region,
        resolution=args.resolution,
        master_path=args.master_path or "maps/",
        save_logs=args.save_logs,
        compute=ComputeConfig(
            days_to_generate=args.days_to_generate,
            criterion=args.criterion,
            max_workers=args.max_workers,
            adaptive_category=(
                default_compute.adaptive_category if args.category_compute is None
                else args.category_compute == "adaptive"
            ),
        ),
        render=RenderConfig(map_mode=args.map_mode),
        conjunction_utc=conjunction,
    )


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.list_methods:
        print(_format_method_listing())
        return 0
    if args.list_criteria:
        print(_format_criteria_listing())
        return 0
    if args.explain_criterion is not None:
        print(_format_criterion_explanation(args.explain_criterion))
        return 0

    if args.command == "render-cache":
        cache = load_visibility_cache_tabular(args.cache_path)
        render_maps_from_cache(
            cache,
            master_path=args.master_path,
            save_logs=bool(args.save_logs),
            perf_report_path=args.perf_report,
        )
        return 0

    try:
        config = _config_from_args(args)
    except ValueError as error:
        parser.error(str(error))
    if args.command in ("cache", "export-csv"):
        export_csv = args.command == "export-csv"
        path = args.csv_path if export_csv else args.cache_path
        output = save_visibility_cache_tabular(build_visibility_cache(config), path)
        print(f"Saved {'CSV export' if export_csv else 'cache'} to {output}")
    else:
        report = args.perf_report
        if args.command == "benchmark" and report is None:
            report = str(config.output_root / "benchmark.json")
        generate_maps(config, perf_report_path=report)
    return 0
