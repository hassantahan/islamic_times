# Workflow Checklist

## Before Coding
1. Identify all affected modules and whether the change crosses Python/C boundaries.
2. Confirm input/output assumptions:
   - Units (degrees, hours, distances)
   - Timezone behavior
   - Dataclass/API contracts
3. Identify any public-facing behavior that could change.

## During Coding
1. Keep changes scoped to the request.
2. Apply strict typing to new/modified Python code.
3. Add or update docstrings for public and non-trivial changed logic (NumPy style for Python).
4. Add comments only where reasoning is not obvious.
5. If C code changes, ensure header declarations and Python integration stay consistent.
6. If C comments/docs are touched, keep wrapper docs and header contracts aligned.

## Validation
1. Run targeted checks for touched paths.
2. For library behavior checks, run `pytest` in `tests/` and relevant script examples (for example `examples/demo.py`) when appropriate.
3. If C-extension code was modified, verify build/import path still works.
4. Confirm no unintended changes outside task scope.
5. Verify venv interpreter type before running build/test commands: `ls .venv/Scripts/python.exe .venv/bin/python 2>/dev/null`.
6. Match shell to venv interpreter type for Python-project commands only (for example `python`, `pip`, `pytest`, `ruff`, extension builds):
   - If `.venv\Scripts\python.exe` exists, use `cmd.exe` for those Python-project commands (including WSL + Windows venv). Example: `cmd.exe /c "cd /d B:\path\to\repo && .venv\Scripts\python.exe -m pytest -q"`.
   - If `.venv/bin/python` exists and `.venv\Scripts\python.exe` does not, use POSIX shell for Python-project commands.
   - Use normal POSIX shell commands for non-Python tasks (for example `rg`, `fd`, `sed`, `git`, file inspection).
7. For broad documentation passes, add/update a tracked audit note summarizing scope, findings, and deferred items.

## Change Summary Template
- What changed:
- Why it changed:
- Public/API impact:
- Edge cases considered:
- Validation performed:
- Remaining risks or follow-up items:

## Release Hygiene
1. For release-bound changes, verify `RELEASE_CHECKLIST.md` items are satisfied.
