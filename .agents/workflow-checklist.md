# Workflow Checklist

## Before Coding

1. Identify all affected modules and whether the change crosses Python/C boundaries.
2. Confirm contracts:
   - Units (degrees, hours, distances)
   - Timezone assumptions
   - Dataclass/API shape
3. Identify any public-facing behavior that may change.

## During Coding

1. Keep changes scoped to the request.
2. Add or update typing, docstrings, and comments for changed non-trivial logic.
3. If C code changes, keep header declarations and Python wrapper contracts aligned.
4. If public behavior changes, update user-facing docs in the same change.

## Validation

1. Run targeted checks for touched paths.
2. Run `pytest` for behavior verification on affected modules.
3. If native code changed, verify build/import path still works.
4. Confirm no unintended edits outside task scope.
5. Select shell for Python-project commands by detected venv type:
   - Check: `ls .venv/Scripts/python.exe .venv/bin/python 2>/dev/null`
   - If `.venv/Scripts/python.exe` exists, use `cmd.exe` for Python-project commands only (`python`, `pip`, `pytest`, `ruff`, extension builds). This includes WSL + Windows venv setups.
   - If `.venv/bin/python` exists and `.venv/Scripts/python.exe` does not, use POSIX shell for Python-project commands.
   - Use POSIX shell for non-Python tasks (`rg`, `git`, `sed`, file inspection).

## Documentation Placement

1. Public docs:
   - `README.md`
   - `CONTRIBUTING.md`
   - `docs/`
2. Internal agent guidance:
   - `.agents/`

## Change Summary Template

- What changed:
- Why it changed:
- Public/API impact:
- Validation performed:
- Remaining risks/follow-up:
