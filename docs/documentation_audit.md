# Documentation Audit - Python and Native C

Date: 2026-02-12  
Scope: `src/islamic_times/`, `src/native/`, `src/native/include/`, `.agents/`

## Summary

This audit focused on documentation quality, consistency, and professional standards
for both Python and native C code. The pass prioritized:

- missing or weak docstrings in public/non-trivial functions,
- ambiguous comments and TODO wording,
- wrapper-level documentation quality for C APIs exposed to Python,
- consistency of project agent guidance for future documentation changes.

No intentional runtime behavior changes were made during this pass.

Second-pass update (same date): additional documentation-only cleanup was applied
to remaining high-visibility Python and C files.

Third-pass update (same date): C-focused contract/comment cleanup was applied to
low-level calculation and datetime files plus related headers.

## Findings and Actions

### 1) Missing/Weak Python Docstrings

Status: **Partially remediated**

Actions completed:

- Added missing docstring to `safe_sun_time` in `src/islamic_times/prayer_times.py`.
- Added missing docstring to `calculate_visibility_shaukat` in `src/islamic_times/moon_equations.py`.
- Upgraded key public `ITLocation` API docstrings in
  `src/islamic_times/islamic_times.py` to cleaner NumPy-style structure.
- Improved key prayer-time function docstrings (`extreme_latitudes`,
  `calculate_prayer_times`) in `src/islamic_times/prayer_times.py`.
- Added short helper docstrings for internal but non-trivial `ITLocation` helpers.
- Normalized public docstrings in `src/islamic_times/time_equations.py` to
  consistent NumPy style and corrected inaccurate return/type descriptions.
- Normalized public docstrings in `src/islamic_times/sun_equations.py`,
  including typo and wording fixes in event-time docs.

### 2) Ambiguous/Non-Professional Comment Wording

Status: **Remediated where identified**

Actions completed:

- Replaced subjective wording such as "personal opinion" with technical intent.
- Replaced "weird issue" wording with explicit day-boundary ordering description.
- Replaced "to be later explained" sentinel notes with explicit C-sentinel behavior.
- Converted a vague Hijri TODO in `src/islamic_times/time_equations.py` into a
  clear limitation note (civil day boundary assumption vs sunset boundary).

### 3) Native C Wrapper/API Documentation Quality

Status: **Improved**

Actions completed:

- Expanded `astro_core` method table doc strings in `src/native/astro_core.c`
  with signature-like guidance and clearer output semantics.
- Improved explanatory comments in `src/native/c_time_equations.c` for units,
  conversion direction, and status behavior.
- Improved `src/native/c_visibilities.c` comments to explain algorithm intent
  and replaced a bare TODO with explicit placeholder status.
- Added concise contract comments in native headers:
  - `src/native/include/c_time_equations.h`
  - `src/native/include/c_sun_equations.h`
  - `src/native/include/c_moon_equations.h`
  - `src/native/include/c_visibilities.h`
- Added/standardized contract and return-code comments in:
  - `src/native/c_sun_equations.c`
  - `src/native/c_moon_equations.c`
- Removed stale wording in native comments (for example placeholder wording
  around loop-term counts).
- C-focused follow-up pass:
  - Normalized section headers and corrected terminology in
    `src/native/c_calculation_equations.c`.
  - Added explicit compatibility note for legacy `compute_equitorial_coordinates`
    naming in both source/header comments.
  - Added function-level contract comments in
    `src/native/include/c_calculation_equations.h`.
  - Added function-level contract comments and clarified comparison semantics in
    `src/native/c_datetime.c`.
  - Added sentinel/struct/function contract comments in
    `src/native/include/c_datetime.h`.
  - Added field-level contract comments for `VisibilityResult` in
    `src/native/include/c_visibilities.h`.

### 4) Agent Guidance Alignment

Status: **Remediated**

Actions completed:

- Updated `.agents/coding-standards.md` to require NumPy-style docstrings for
  new/modified Python docs and clearer C contract-comment expectations.
- Updated `.agents/workflow-checklist.md` to include documentation-specific
  workflow checkpoints for broad doc passes.
- Clarified shell policy for WSL + Windows venv workflows: use `cmd.exe` for
  Python-project commands only, and POSIX shell for non-Python tasks.

## Deferred / Follow-up Items

The following are quality opportunities identified during review and left for
subsequent focused cleanup passes:

1. Normalize remaining legacy Python docstrings in:
   - `src/islamic_times/moon_equations.py`
2. Reduce high-volume `# type: ignore` usage in
   `src/islamic_times/prayer_times.py` with tighter typing refinements.
3. Continue harmonizing native C comment style in remaining files not heavily
   touched in this pass (optional polish; no high-priority gaps identified in
   touched C files after this follow-up).

## Validation Notes

- This pass was documentation/comment focused and did not intentionally modify
  algorithm behavior.
- Runtime validation should still be executed after large doc-only refactors to
  ensure no accidental logic edits were introduced.

### Latest Validation Snapshot (2026-02-12)

- `cmd.exe /c ".venv\Scripts\python.exe -m py_compile src\islamic_times\sun_equations.py src\islamic_times\time_equations.py"` passed.
- `cmd.exe /c ".venv\Scripts\python.exe -m pytest -q"` passed (full test suite).
- `ruff check` on the touched Python files still reports pre-existing baseline
  style issues; not remediated in this documentation-focused pass.
- `cmd.exe /c ".venv\Scripts\python.exe setup.py build_ext --inplace"` passed
  after C-comment/header updates.
- `cmd.exe /c ".venv\Scripts\python.exe -m pytest -q tests\test_native_validation.py tests\test_transit_day_normalization.py tests\test_calculation_equations.py"`
  passed (native-focused validation subset).

## Completion Status

Status: **Audit complete for current scope**

- Python and C documentation/docstring/comment baselines were reviewed.
- High-priority C documentation-contract gaps identified during this review are
  now remediated in the touched files.
- Remaining follow-up items are optional quality polish, not blockers.

### Issue 3 Verification Record (Resolved)

Issue 3 corresponds to the B-series loop bound in native moon nutation logic.
The active source now uses:

- `src/native/c_moon_equations.c`: `for (size_t i = 0; i < n_args; i++)`

for the B-series iteration in `moon_nutation`.

Validation sample (same reproducible input set):

- Inputs: `jde=2460680.5`, `delta_t=69.5`, `lat=40.7128`, `lon=-74.0060`,
  `elev_m=10.0`, `temp_c=15.0`, `pressure_kpa=101.0`.
- Post-fix native vs Python deltas:
  - `delta_sum_l=-0.000000002328`
  - `delta_sum_b=0.000000000349`
  - `delta_sum_r=0.000000007451`

These are floating-point rounding-level differences and confirm parity for the
previously failing B-series term summation path.
