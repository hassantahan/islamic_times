# Contributing

## 1. Development Setup

1. Create and activate a Python 3.10+ virtual environment.
2. Install package + development tools:

```bash
pip install -e ".[dev]"
```

3. If working on mapper code, install optional mapper dependencies:

```bash
pip install -e ".[map]"
```

## 2. Shell and Interpreter Policy

Use POSIX shell for non-Python commands (`rg`, `git`, `sed`, file inspection).
For Python-project commands only (`python`, `pip`, `pytest`, `ruff`, extension builds), pick the shell by venv type:

1. Verify venv interpreter type:

```bash
ls .venv/Scripts/python.exe .venv/bin/python 2>/dev/null
```

2. If `.venv/Scripts/python.exe` exists, run Python-project commands with `cmd.exe`.
This includes WSL sessions using a Windows venv.

```bat
cmd.exe /c "cd /d B:\path\to\repo && .venv\Scripts\python.exe -m pytest -q"
```

3. If `.venv/bin/python` exists and `.venv/Scripts/python.exe` does not, use normal POSIX Python commands.

## 3. Validation Before PR

Run checks relevant to your changes. Typical baseline:

```bash
ruff check src tests
pytest \
  --cov=islamic_times.islamic_times \
  --cov=islamic_times.it_dataclasses \
  --cov=islamic_times.calculation_equations \
  --cov=islamic_times.time_equations \
  --cov-report=term-missing \
  --cov-fail-under=85
```

If native C files changed, also verify extension build/import path:

```bash
python setup.py build_ext --inplace
```

## 4. Pull Request Expectations

- Keep each PR scoped to one primary topic.
- Describe behavior impact and compatibility impact.
- Include exact validation commands executed and outcomes.
- Update public docs (`README.md`, `CONTRIBUTING.md`, `docs/`) when user-facing behavior changes.

## 5. Release Workflow

For release-bound changes, follow `RELEASE_CHECKLIST.md`.
