# Contributing

## Development Setup
1. Create and activate a Python 3.10+ virtual environment.
2. Install project and development dependencies:

```bash
pip install -e ".[dev]"
```

3. If you work on map generation tools, install optional mapping dependencies:

```bash
pip install -e ".[map]"
```

### Windows examples

PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
```

cmd.exe:

```bat
.venv\Scripts\activate.bat
python -m pip install -e ".[dev]"
```

## Running Checks
Run these checks before opening a pull request:

```bash
ruff check tests
pytest \
  --cov=islamic_times.islamic_times \
  --cov=islamic_times.it_dataclasses \
  --cov=islamic_times.calculation_equations \
  --cov=islamic_times.time_equations \
  --cov-report=term-missing \
  --cov-fail-under=85
```

## C Extension Notes
- C extension sources are in `src/native/*.c`.
- Headers are in `src/native/include/*.h`.
- Keep function signatures and units synchronized between C and Python call sites.

## Pull Requests
- Keep changes scoped to one topic.
- Include a short summary of behavior impact.
- Mention validation steps and results.

## Releases
See `RELEASE_CHECKLIST.md` for tag-and-publish steps.
