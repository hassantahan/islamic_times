"""CLI for mapper generation and benchmarking."""

from __future__ import annotations

import argparse
from datetime import datetime

from .config import ComputeConfig, MapperConfig, RenderConfig
from .pipeline import generate_maps
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
    for method in pt.DEFAULT_PRAYER_METHODS:
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
    parser = MapperArgumentParser(description="Generate new-moon visibility maps")
    subparsers = parser.add_subparsers(dest="command", required=False)

    def add_common_args(target: argparse.ArgumentParser) -> None:
        target.add_argument("--date", type=str, default=None, help="ISO datetime for month/year (e.g. 2025-01-01T00:00:00)")
        target.add_argument("--master_path", type=str, default="maps/", help="Directory where maps are written")
        target.add_argument("--total_months", type=int, default=1)
        target.add_argument("--map_region", type=str, default="WORLD")
        target.add_argument("--map_mode", type=str, default="category", choices=("raw", "category"))
        target.add_argument("--resolution", type=int, default=300)
        target.add_argument("--days_to_generate", type=int, default=3)
        target.add_argument("--criterion", type=int, default=1, choices=(0, 1, 2))
        target.add_argument("--save_logs", action="store_true")
        target.add_argument("--max_workers", type=int, default=None)
        target.add_argument("--perf_report", type=str, default=None, help="Optional JSON output path for timing report")

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
            default=None,
            help="Print details for one visibility criterion and exit.",
        )

    common = argparse.ArgumentParser(add_help=False)
    add_common_args(common)
    add_common_args(parser)
    add_info_args(common)
    add_info_args(parser)

    subparsers.add_parser("generate", parents=[common], help="Generate mapper outputs")
    subparsers.add_parser("benchmark", parents=[common], help="Generate outputs and emit benchmark JSON")
    return parser


def _config_from_args(args: argparse.Namespace) -> MapperConfig:
    if args.criterion not in _CRITERION_INFO:
        raise ValueError(
            f"Invalid criterion '{args.criterion}'. Supported values are 0 (Odeh), 1 (Yallop), and 2 (Shaukat)."
        )
    date = datetime.fromisoformat(args.date) if args.date else datetime.now()
    return MapperConfig(
        date=date,
        total_months=args.total_months,
        map_region=args.map_region,
        resolution=args.resolution,
        master_path=args.master_path,
        save_logs=args.save_logs,
        compute=ComputeConfig(
            days_to_generate=args.days_to_generate,
            criterion=args.criterion,
            max_workers=args.max_workers,
        ),
        render=RenderConfig(map_mode=args.map_mode),
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

    command = args.command or "generate"

    cfg = _config_from_args(args)
    perf_path = args.perf_report if command in ("benchmark", "generate") else None
    generate_maps(cfg, perf_report_path=perf_path)
    return 0
