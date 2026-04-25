from __future__ import annotations

from datetime import datetime, timezone

import islamic_times.astro_core as fast_astro
import numpy as np
import pytest

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import PUBLIC_SCHEMA_VERSION, BatchVisibilities


def test_batch_visibilities_raw_returns_typed_matrix() -> None:
    result = ITLocation.batch_visibilities(
        latitudes=[43.651070, 40.7128],
        longitudes=[-79.347015, -74.0060],
        date=datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),
        days=2,
        criterion=1,
        output="raw",
    )

    assert isinstance(result, BatchVisibilities)
    assert result.criterion == "Yallop"
    assert result.shape == (2, 2)
    assert isinstance(result.values[0][0], float)


def test_batch_visibilities_classification_and_code_modes() -> None:
    date = datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc)
    labels = ITLocation.batch_visibilities(
        latitudes=[43.651070],
        longitudes=[-79.347015],
        date=date,
        days=1,
        criterion=2,
        output="classification",
    )
    codes = ITLocation.batch_visibilities(
        latitudes=[43.651070, 40.7128],
        longitudes=[-79.347015, -74.0060],
        date=date,
        days=1,
        criterion=2,
        output="code",
    )

    assert labels.output == "classification"
    assert labels.criterion == "Shaukat"
    assert isinstance(labels.values[0][0], str)

    assert codes.output == "code"
    assert codes.shape == (2, 1)
    assert isinstance(codes.values[0][0], int)
    assert max(codes.values[0][0], codes.values[1][0]) <= 9


def test_batch_visibilities_matches_native_batch_result_for_raw_mode() -> None:
    date = datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc)
    batch = ITLocation.batch_visibilities(
        latitudes=[43.651070],
        longitudes=[-79.347015],
        date=date,
        days=2,
        criterion=1,
        utc_offset=0.0,
        elevation=10.0,
        temperature=15.0,
        pressure=101.325,
        output="raw",
    )
    native = fast_astro.compute_visibilities_batch(
        np.array([43.651070], dtype=np.float64),
        np.array([-79.347015], dtype=np.float64),
        date,
        2,
        1,
        0.0,
        10.0,
        15.0,
        101.325,
        "r",
    )

    assert batch.shape == (1, 2)
    assert batch.values[0][0] == pytest.approx(float(native[0]), rel=0.0, abs=1e-12)
    assert batch.values[0][1] == pytest.approx(float(native[1]), rel=0.0, abs=1e-12)


def test_batch_visibilities_validates_coordinate_lengths_and_output() -> None:
    date = datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc)

    with pytest.raises(ValueError, match="same length"):
        ITLocation.batch_visibilities(
            latitudes=[43.651070, 40.7128],
            longitudes=[-79.347015],
            date=date,
            output="raw",
        )

    with pytest.raises(ValueError, match="must be one of"):
        ITLocation.batch_visibilities(
            latitudes=[43.651070],
            longitudes=[-79.347015],
            date=date,
            output="labels",  # type: ignore[arg-type]
        )


def test_batch_visibilities_to_dict_contract() -> None:
    result = ITLocation.batch_visibilities(
        latitudes=[43.651070],
        longitudes=[-79.347015],
        date=datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),
        days=1,
        criterion=1,
        output="raw",
    )
    payload = result.to_dict()

    assert payload["schema_version"] == PUBLIC_SCHEMA_VERSION
    assert payload["type"] == "BatchVisibilities"
    assert payload["criterion"] == "Yallop"
    assert payload["output"] == "raw"
    assert payload["shape"] == [1, 1]
