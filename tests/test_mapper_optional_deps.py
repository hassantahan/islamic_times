from __future__ import annotations

import pytest

import mapper


def test_mapper_dependency_guard() -> None:
    if mapper._MAPPER_IMPORT_ERROR is None:
        mapper.require_mapper_dependencies()
    else:
        with pytest.raises(ImportError):
            mapper.require_mapper_dependencies()
