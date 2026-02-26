"""Category label and color palettes for crescent visibility maps."""

from __future__ import annotations

from typing import Dict, List, Tuple

CategoryEntry = Tuple[str, str, float]

_SPECIAL_PREFIX: list[CategoryEntry] = [
    ("Moonset before the new moon.", "#141414", 1.0),
    ("Moonset before sunset.", "#393a3c", 1.0),
    ("Moonset & Sunset don't exist.", "#233342", 1.0),
    ("Sunset doesn't exist.", "#35526b", 1.0),
    ("Moonset doesn't exist.", "#5b6875", 1.0),
]

_CRITERION_LABELS: dict[int, list[CategoryEntry]] = {
    0: _SPECIAL_PREFIX
    + [
        ("D: Crescent is not visible even by optical aid.", "#807f80", 0.15),
        ("C: Crescent is visible by optical aid only.", "#B89D18", 1.0),
        ("B: Crescent is visible by optical aid, and it could be seen by naked eyes.", "#74b818", 1.0),
        ("A: Crescent is visible by naked eyes.", "#1BB818", 1.0),
    ],
    1: _SPECIAL_PREFIX
    + [
        ("F: Not visible; below the Danjon limit.", "#807f80", 0.15),
        ("E: Not visible with a [conventional] telescope.", "#B81818", 1.0),
        ("D: Will need optical aid to find crescent.", "#e3d61b", 1.0),
        ("C: May need optical aid to find crescent.", "#89d518", 1.0),
        ("B: Visible under perfect conditions.", "#54b818", 1.0),
        ("A: Easily visible.", "#1bdf18", 1.0),
    ],
}


def category_entries(criterion: int) -> list[CategoryEntry]:
    """Return ordered category entries for a supported criterion."""
    if criterion not in _CRITERION_LABELS:
        raise ValueError(f"Unsupported criterion '{criterion}'.")
    return list(_CRITERION_LABELS[criterion])


def category_labels(criterion: int) -> List[str]:
    """Return ordered category labels for palette and integer-code mapping."""
    return [label for label, _, _ in category_entries(criterion)]


def category_code_map(criterion: int) -> Dict[str, int]:
    """Map classification strings to compact integer codes."""
    return {label: idx for idx, label in enumerate(category_labels(criterion))}

