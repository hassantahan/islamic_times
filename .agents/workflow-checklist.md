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
3. Add or update docstrings for public and non-trivial changed logic.
4. Add comments only where reasoning is not obvious.
5. If C code changes, ensure header declarations and Python integration stay consistent.

## Validation
1. Run targeted checks for touched paths.
2. For library behavior checks, run `pytest` in `tests/` and relevant script examples (for example `examples/demo.py`) when appropriate.
3. If C-extension code was modified, verify build/import path still works.
4. Confirm no unintended changes outside task scope.

## Change Summary Template
- What changed:
- Why it changed:
- Public/API impact:
- Edge cases considered:
- Validation performed:
- Remaining risks or follow-up items:

## Release Hygiene
1. For release-bound changes, verify `RELEASE_CHECKLIST.md` items are satisfied.
