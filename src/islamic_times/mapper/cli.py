"""CLI for mapper generation and benchmarking."""

from __future__ import annotations

import argparse
from datetime import datetime

from .config import ComputeConfig, MapperConfig, RenderConfig
from .pipeline import generate_maps


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate new-moon visibility maps")
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

    common = argparse.ArgumentParser(add_help=False)
    add_common_args(common)
    add_common_args(parser)

    subparsers.add_parser("generate", parents=[common], help="Generate mapper outputs")
    subparsers.add_parser("benchmark", parents=[common], help="Generate outputs and emit benchmark JSON")
    return parser


def _config_from_args(args: argparse.Namespace) -> MapperConfig:
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
    command = args.command or "generate"

    cfg = _config_from_args(args)
    perf_path = args.perf_report if command in ("benchmark", "generate") else None
    generate_maps(cfg, perf_report_path=perf_path)
    return 0
