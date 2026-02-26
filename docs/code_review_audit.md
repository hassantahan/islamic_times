# Code Review Audit - Refactor, Restructure, and Performance

Date: 2026-02-12  
Scope: `src/islamic_times/`, `src/native/`, `src/native/include/`, `tests/`  
Out of scope in this pass: `mapper.py` and mapping workflow internals.

## Summary

This audit reviewed Python orchestration, native C kernels/wrappers, and tests with three goals:

1. Identify correctness risks that can cause wrong outputs or long-running instability.
2. Identify maintainability/organization gaps that slow future development.
3. Identify practical performance opportunities with measurable impact.

Overall assessment:

- Numerical core performance is strong, but there are critical correctness and API-contract gaps around native wrappers.
- Module organization has drifted into mixed responsibilities (active API + deprecated legacy implementations in same modules).
- Type-safety and error-model consistency are weak in prayer/extreme-latitude logic, creating broad exception handling and many `# type: ignore` paths.

## Method and Baseline Evidence

### Static review pass

- Reviewed full core Python modules:
  - `src/islamic_times/islamic_times.py`
  - `src/islamic_times/prayer_times.py`
  - `src/islamic_times/sun_equations.py`
  - `src/islamic_times/moon_equations.py`
  - `src/islamic_times/time_equations.py`
  - `src/islamic_times/calculation_equations.py`
  - `src/islamic_times/it_dataclasses.py`
- Reviewed full native C stack and headers:
  - `src/native/*.c`
  - `src/native/include/*.h`
- Reviewed tests in `tests/` for coverage and regression depth.

### Runtime probes (Windows venv via `cmd.exe`)

- `ITLocation.__init__(find_local_tz=False, auto_calculate=True)`: ~0.124 ms/op
- `ITLocation.calculate_prayer_times()`: ~0.029 ms/op
- `ITLocation.visibilities(days=3)`: ~0.069 ms/op
- `astro_core.compute_visibilities_batch` with 200 coordinates: ~4.792 ms/op
- `ITLocation.__init__(find_local_tz=True)`: ~249.280 ms/op

Interpretation: local timezone lookup dominates runtime when enabled.

### Import-time probe

`python -X importtime tmp_import_probe.py` showed:

- `islamic_times.islamic_times` cumulative import: ~169 ms
- `islamic_times.time_equations` import path dominates and eagerly loads heavy timezonefinder/h3 stack.

### Targeted correctness/performance probes

- Native API `criterion=2` (Shaukat) currently returns constant `0.0` q-values and repeated class `"B: Visible under perfect conditions."`.
- Memory probe indicates growth when calling `astro_core.compute_moon` in loops; growth persists after `gc.collect()`.
- Scalar helper benchmark: `calculation_equations.sin` is ~13.9x slower than `math.sin(math.radians(x))` for scalar loops.

## Prioritized Backlog

Legend:

- Severity: `High`, `Medium`, `Low`
- Effort: `S` (small), `M` (medium), `L` (large)
- Risk: implementation risk, not business priority

---

### CR-001 - Native memory leak in `compute_moon` wrapper

- Severity: `High`
- Category: `correctness`, `performance`, `stability`
- Evidence:
  - `src/native/c_moon_equations.c:669`
  - `src/native/c_moon_equations.c:671`
  - `src/native/c_moon_equations.c:672`
  - `src/native/c_moon_equations.c:673`
- Why it matters:
  - Long-running processes repeatedly calling lunar computations can steadily consume memory.
  - Probe showed monotonic growth after repeated calls and GC.
- Root cause (specific):
  - `PyTuple_Pack` is used with inline `PyFloat_FromDouble(...)` temporaries for `sum_l/sum_b/sum_r`, creating unbalanced ref handling.
- Recommended change:
  - Build tuple explicitly with `PyTuple_New(4)` and `PyTuple_SET_ITEM`, or create temporary refs and `Py_DECREF` after `PyTuple_Pack`.
  - Add a dedicated native wrapper leak regression test.
- Expected impact:
  - Prevent process memory growth under repeated moon computations.
- Effort: `S`
- Risk: `Low`
- Validation:
  - Repeat the current memory probe and confirm flat traced memory after warmup.
  - Add a stress test loop in CI (smaller iteration count, leak-threshold assertion).

---

### CR-002 - Shaukat criterion is exposed as valid but not implemented

- Severity: `High`
- Category: `correctness`, `api-clarity`
- Evidence:
  - Placeholder implementation: `src/native/c_visibilities.c:29`
  - Criterion accepted by wrappers: `src/native/c_visibilities.c:223`, `src/native/c_visibilities.c:343`
  - Native output labels criterion as `"Shaukat"`: `src/native/c_visibilities.c:286`
  - Public Python facade rejects `2`: `src/islamic_times/islamic_times.py:652`
- Why it matters:
  - Direct native API callers receive plausible-looking but invalid results.
  - Contract inconsistency between public facade and native extension can cause integration errors.
- Recommended change:
  - Short term: reject criterion `2` in all native wrappers with explicit `NotImplementedError` semantics.
  - Medium term: implement actual Shaukat model and add parity/expected-value tests.
- Expected impact:
  - Eliminates silent false outputs from native API.
- Effort: `S` (reject path) / `M` (full implementation)
- Risk: `Low` (reject path) / `Medium` (full implementation)
- Validation:
  - Add tests asserting criterion `2` currently raises until implemented.
  - When implemented, add known-case dataset tests.

---

### CR-003 - Broad exception translation masks real defects

- Severity: `High`
- Category: `correctness`, `maintainability`
- Evidence:
  - `src/islamic_times/sun_equations.py:528`
  - `src/islamic_times/moon_equations.py:642`
  - `src/islamic_times/prayer_times.py:367`
  - `src/islamic_times/prayer_times.py:396`
  - `src/islamic_times/prayer_times.py:499`
- Why it matters:
  - `except Exception` converts programming errors and data-contract bugs into domain `ArithmeticError`/`inf`, obscuring root causes.
  - Debugging and incident response become significantly harder.
- Recommended change:
  - Narrow exception clauses to expected domain failures only.
  - Preserve original exception types for programmer/configuration errors.
  - Replace sentinel fallback in non-domain failures with explicit error propagation.
- Expected impact:
  - Faster defect localization; fewer silent wrong outputs.
- Effort: `M`
- Risk: `Medium` (may expose currently hidden failures)
- Validation:
  - Add tests verifying invalid types/ranges raise specific exception classes.
  - Add negative tests to ensure domain non-existence still maps to intended error paths.

---

### CR-004 - Sentinel-heavy time model causes fragile control flow

- Severity: `High`
- Category: `architecture`, `maintainability`
- Evidence:
  - Mixed `datetime | str | float` prayer value type: `src/islamic_times/it_dataclasses.py:413`
  - Sentinel checks and string parsing in core logic: `src/islamic_times/prayer_times.py:229`, `src/islamic_times/prayer_times.py:516`
  - `time_midpoint` sentinel checks: `src/islamic_times/time_equations.py:452`
- Why it matters:
  - Business logic depends on runtime union disambiguation and string content checks.
  - Drives high `# type: ignore` count and broad exception blocks.
- Recommended change:
  - Introduce structured event result type, e.g.:
    - `PrayerValue(value: datetime | None, status: PrayerStatus, message: str | None)`
  - Keep public `Prayer.time` compatibility layer temporarily, but internal computations should use typed status objects.
- Expected impact:
  - Lower complexity in `extreme_latitudes`, fewer type-ignore paths, safer refactors.
- Effort: `L`
- Risk: `Medium` (touches many call sites)
- Validation:
  - Add internal invariants tests (no string parsing in algorithm layer).
  - Snapshot tests for public string formatting backward compatibility.

---

### CR-005 - Deprecated Python implementations are co-located with active runtime paths

- Severity: `Medium`
- Category: `organization`, `performance`
- Evidence:
  - Deprecated heavy routines remain in active modules:
    - `src/islamic_times/sun_equations.py:214`
    - `src/islamic_times/moon_equations.py:323`
    - `src/islamic_times/time_equations.py:76`
  - Import-time weight dominated by modules eagerly imported in normal usage.
- Why it matters:
  - Increases startup/import latency and maintenance surface.
  - Blurs which implementation is canonical.
- Recommended change:
  - Move deprecated pure-Python algorithms to `islamic_times._legacy_py_impl`.
  - Keep thin compatibility wrappers that lazy-import legacy paths only when called.
  - Mark deprecation timeline in docs/changelog.
- Expected impact:
  - Cleaner architecture and lower default import/runtime overhead.
- Effort: `M`
- Risk: `Low`
- Validation:
  - Compare `-X importtime` before/after.
  - Run compatibility tests for legacy wrappers.

---

### CR-006 - Timezone lookup path is materially expensive and un-cached

- Severity: `Medium`
- Category: `performance`
- Evidence:
  - Heavy object construction every call: `src/islamic_times/time_equations.py:409`
  - Eager import of timezonefinder stack: `src/islamic_times/time_equations.py:19`
  - Measured `find_local_tz=True` init cost ~249 ms/op.
- Why it matters:
  - Dominates latency for users who enable local timezone detection.
- Recommended change:
  - Lazy import timezonefinder and `pytz` only in `find_utc_offset`.
  - Reuse a module-level singleton `TimezoneFinder`.
  - Add LRU cache keyed by rounded `(lat, lon, date)` where acceptable.
- Expected impact:
  - Significant latency reduction for repeated timezone resolution.
- Effort: `S/M`
- Risk: `Low`
- Validation:
  - Add microbenchmark guard for repeated calls.
  - Add correctness tests around DST boundaries.

---

### CR-007 - Scalar trig helpers use NumPy and are slow in scalar paths

- Severity: `Medium`
- Category: `performance`
- Evidence:
  - `src/islamic_times/calculation_equations.py:14`
  - `src/islamic_times/calculation_equations.py:24`
  - Measured scalar loop slowdown ~13.9x for `ce.sin` vs `math.sin(math.radians(x))`.
- Why it matters:
  - For scalar-heavy Python fallback paths, this introduces avoidable overhead.
- Recommended change:
  - Use `math.sin/math.cos` for scalar helpers.
  - Keep vectorized NumPy variants only where arrays are passed.
- Expected impact:
  - Lower overhead in scalar helper call sites.
- Effort: `S`
- Risk: `Low`
- Validation:
  - Keep numeric parity tests for helper outputs.
  - Benchmark scalar helper loops pre/post.

---

### CR-008 - `extreme_latitudes` is monolithic and difficult to reason about

- Severity: `Medium`
- Category: `architecture`, `maintainability`
- Evidence:
  - `src/islamic_times/prayer_times.py:170` through `src/islamic_times/prayer_times.py:421`
  - Extensive type-ignore/branching and mixed responsibilities.
- Why it matters:
  - High regression risk for any change in extreme-latitude behavior.
  - Hard to isolate and test sub-rules.
- Recommended change:
  - Split into strategy functions per rule (`NEARESTLAT`, `MIDDLENIGHT`, `ONESEVENTH`, `ANGLEBASED`).
  - Introduce intermediate typed context object (precomputed events/angles).
  - Make each strategy pure and unit-testable.
- Expected impact:
  - Reduced complexity and safer rule-specific changes.
- Effort: `L`
- Risk: `Medium`
- Validation:
  - Golden tests per strategy and latitude band.
  - Existing integration tests must remain green.

---

### CR-009 - Duplicated native solver logic across sun/moon paths

- Severity: `Medium`
- Category: `organization`, `maintainability`, `performance`
- Evidence:
  - Solar duplication:
    - `src/native/c_sun_equations.c:596`
    - `src/native/c_sun_equations.c:737`
  - Lunar duplication:
    - `src/native/c_moon_equations.c:855`
    - `src/native/c_moon_equations.c:939`
- Why it matters:
  - Bug fixes need to be mirrored in multiple places.
  - Increases risk of drift and inconsistent behavior.
- Recommended change:
  - Extract shared internal interpolation/event iteration helpers.
  - Parameterize with body-specific coordinate callbacks.
- Expected impact:
  - Lower maintenance cost and fewer divergence bugs.
- Effort: `M/L`
- Risk: `Medium`
- Validation:
  - Re-run transit/day-normalization tests.
  - Add parity tests between old/new paths during migration.

---

### CR-010 - Unbounded `while(1)` day-alignment loops in native code

- Severity: `Medium`
- Category: `correctness`, `maintainability`
- Evidence:
  - `src/native/c_sun_equations.c:687`
  - `src/native/c_sun_equations.c:831`
  - `src/native/c_moon_equations.c:955`
- Why it matters:
  - Infinite-loop risk is low but non-zero under unexpected edge conditions.
  - No explicit iteration cap or diagnostic if convergence fails.
- Recommended change:
  - Add bounded iteration cap and explicit failure status codes.
  - Return diagnostic error for non-convergence.
- Expected impact:
  - Safer failure behavior and easier debugging.
- Effort: `S/M`
- Risk: `Low`
- Validation:
  - Add tests that assert controlled failure on forced non-convergent fixtures.

---

### CR-011 - API inconsistency between facade and native criterion support

- Severity: `Medium`
- Category: `api-clarity`
- Evidence:
  - Facade only permits `criterion in (0, 1)`: `src/islamic_times/islamic_times.py:652`
  - Native wrappers permit `0..2`: `src/native/c_visibilities.c:223`, `src/native/c_visibilities.c:343`
- Why it matters:
  - Users of native API and high-level API observe inconsistent contracts.
- Recommended change:
  - Define single source of truth for criterion enum and enforce consistently.
  - Export a shared constant map at Python layer.
- Expected impact:
  - Fewer integration surprises.
- Effort: `S`
- Risk: `Low`
- Validation:
  - Cross-layer API contract tests for accepted/rejected values.

---

### CR-012 - Global native type cache uses process-exit cleanup pattern

- Severity: `Medium`
- Category: `professional-standards`, `maintainability`
- Evidence:
  - Global type pointers: `src/native/astro_core.c:9`
  - `Py_AtExit` cleanup registration: `src/native/astro_core.c:138`
- Why it matters:
  - Global state patterns are less robust in multi-interpreter contexts.
  - Harder to evolve toward modern module-state patterns.
- Recommended change:
  - Migrate to per-module state (`PyModule_GetState`) and avoid process-global caches.
- Expected impact:
  - Cleaner lifecycle management and future-proofing.
- Effort: `M`
- Risk: `Medium`
- Validation:
  - Module init/teardown tests and repeated import/unload stress.

---

### CR-013 - Sequence parsing in native moon wrappers is overly strict and partially unsafe

- Severity: `Medium`
- Category: `correctness`, `api-robustness`
- Evidence:
  - `src/native/c_moon_equations.c:812`
  - `src/native/c_moon_equations.c:815`
  - `src/native/c_moon_equations.c:1009`
  - `src/native/c_moon_equations.c:1012`
- Why it matters:
  - Rejects integer values though they are numerically valid.
  - Does not first guard `PySequence_GetItem` return values before type checks.
- Recommended change:
  - Use `PyFloat_AsDouble` with explicit error checks (accept `int` and `float`).
  - Guard null returns before type/value checks.
- Expected impact:
  - Better API ergonomics and safer error handling.
- Effort: `S`
- Risk: `Low`
- Validation:
  - Add tests for mixed numeric input sequences and malformed sequences.

---

### CR-014 - `find_utc_offset` does not handle missing timezone resolution explicitly

- Severity: `Medium`
- Category: `correctness`
- Evidence:
  - `src/islamic_times/time_equations.py:412`
  - `src/islamic_times/time_equations.py:415`
- Why it matters:
  - If timezone lookup returns `None`, downstream failure messages are unclear.
- Recommended change:
  - Validate `timezone_str` and raise explicit `ValueError` with coordinates/date context.
  - Add fallback strategy (optional): nearest valid zone or UTC with warning.
- Expected impact:
  - Better error quality and operational clarity.
- Effort: `S`
- Risk: `Low`
- Validation:
  - Add tests with forced invalid/unresolvable coordinate scenarios.

---

### CR-015 - Repeated linear method lookup in `set_prayer_method`

- Severity: `Low`
- Category: `performance`, `organization`
- Evidence:
  - `src/islamic_times/islamic_times.py:384`
- Why it matters:
  - Minor overhead, but unnecessary repeated scans.
- Recommended change:
  - Build one-time alias map `{normalized_key: PrayerMethod}` at module import.
- Expected impact:
  - Small speedup and cleaner method selection code.
- Effort: `S`
- Risk: `Low`
- Validation:
  - Existing method-selection tests should remain unchanged.

---

### CR-016 - Repetitive stale-state guards in `ITLocation` accessors

- Severity: `Low`
- Category: `organization`
- Evidence:
  - `src/islamic_times/islamic_times.py:545`
  - `src/islamic_times/islamic_times.py:562`
  - `src/islamic_times/islamic_times.py:592`
  - `src/islamic_times/islamic_times.py:609`
- Why it matters:
  - Duplicates policy logic and message variants.
- Recommended change:
  - Extract a single `_ensure_state_ready(target: Literal[...])` guard helper.
- Expected impact:
  - Lower duplication and easier policy changes.
- Effort: `S`
- Risk: `Low`
- Validation:
  - State-transition tests remain green.

---

### CR-017 - Sun core receives weather inputs that are effectively unused

- Severity: `Low`
- Category: `api-clarity`
- Evidence:
  - `compute_sun_result` signature includes `temperature/pressure`: `src/native/c_sun_equations.c:253`
  - Apparent altitude currently equals true altitude: `src/native/c_sun_equations.c:364`
- Why it matters:
  - API suggests atmospheric dependency not currently applied.
- Recommended change:
  - Either re-enable physically justified refraction path with tests, or simplify signature/docs to avoid misleading parameters.
- Expected impact:
  - Better contract clarity.
- Effort: `S/M`
- Risk: `Low`
- Validation:
  - API docs and wrapper tests updated accordingly.

---

### CR-018 - Legacy sentinel constants reduce readability and auditability

- Severity: `Low`
- Category: `professional-standards`
- Evidence:
  - `src/native/include/c_moon_equations.h:7`
  - `src/islamic_times/moon_equations.py:475`
- Why it matters:
  - Magic values (`-123456.0`) are fragile and opaque.
- Recommended change:
  - Replace with explicit enum/status flags or dedicated structs for optional nutation reuse.
- Expected impact:
  - Cleaner contracts and fewer accidental sentinel collisions.
- Effort: `M`
- Risk: `Low`
- Validation:
  - Wrapper parsing tests and transit/moontime regression suite.

## Public API / Interface Impact (Proposed)

No changes were implemented in this review pass. Proposed API-impacting changes for later implementation include:

1. Criterion handling:
   - Align native and Python contracts for visibility criterion support.
2. Prayer/event typing:
   - Move from sentinel unions toward structured status objects (with compatibility layer).
3. Legacy function placement:
   - Relocate deprecated pure-Python algorithms behind lazy compatibility wrappers.

Each of these should ship behind a compatibility plan with explicit versioned deprecation notes.

## Test Gaps and Recommended Additions

### Correctness regressions

1. Native memory stability test:
   - Repeated `compute_moon` calls should not produce sustained traced-memory growth after warmup.
2. Criterion contract test:
   - Criterion `2` behavior must be explicit (rejected until implemented, or verified if implemented).
3. Error-class test matrix:
   - Ensure domain non-existence errors are distinct from programmer/configuration errors.

### Behavioral edge cases

1. Extreme-latitude strategy matrix:
   - Parameterized tests by rule and latitude band.
2. Timezone resolution failure paths:
   - Explicit tests for unresolved coordinates/invalid timezone outcomes.

### Performance guardrails

1. Optional benchmark smoke tests:
   - `find_utc_offset` repeated calls with/without cache.
   - Scalar helper microbench (`ce.sin`/`ce.cos`) if retained.
2. Import-time budget check:
   - Track startup cost trend when touching `time_equations` and top-level imports.

## Implementation Roadmap

### Phase 1 (Immediate: correctness and safety)

1. Fix `compute_moon` wrapper leak (`CR-001`).
2. Disable or explicitly reject unimplemented Shaukat path (`CR-002`, `CR-011`).
3. Narrow broad exception translation (`CR-003`).

### Phase 2 (Near-term: performance + maintainability)

1. Lazy import and cache timezone lookup dependencies (`CR-006`, `CR-014`).
2. Replace scalar NumPy trig helpers with scalar math implementations (`CR-007`).
3. Tighten sequence parsing robustness in moon wrappers (`CR-013`).
4. Add bounded iteration protections for native day-alignment loops (`CR-010`).

### Phase 3 (Structural refactor)

1. Split `extreme_latitudes` into strategy functions and typed contexts (`CR-004`, `CR-008`).
2. Consolidate duplicated native event solver logic (`CR-009`).
3. Separate legacy Python algorithms from active runtime modules (`CR-005`).
4. Consider module-state migration in native extension (`CR-012`).

Phase 3 execution status (implemented):

1. Done: `extreme_latitudes` was split into strategy-oriented helpers with typed intermediate structures.
2. Done: shared native event-solver helpers now back solar/lunar transit and rise/set search flows to remove duplicated iteration logic.
3. Done: legacy Python implementations are now isolated under `src/islamic_times/_legacy_py_impl/` and loaded lazily by deprecated facades.
4. Done: native extension type caches now live in per-module state (`m_size` + `traverse`/`clear`/`free`) and no longer depend on process-exit cleanup.

### Phase 4 (Stabilization + optimization)

Phase 4 execution status (implemented):

1. Added targeted native hot-loop optimizations in shared event-solver helpers (`src/native/c_event_solver.c`):
   - Precompute invariant trigonometric values in altitude refinement loops.
   - Remove redundant degree-radian reconversion in `delta_m` updates.
   - Use direct day-level datetime comparisons in day-search helper to reduce repeated `day_of_year` branching.
2. Reduced Python C-API call overhead in native wrappers:
   - Use `PyObject_CallOneArg` where available for one-argument dataclass constructors in sun/moon wrappers.
3. Reduced timezone lookup overhead in Python time helpers (`src/islamic_times/time_equations.py`):
   - Added cached `pytz.timezone` factory and per-timezone object cache.
4. Added reproducible performance probe script:
   - `examples/perf_phase4_probe.py`

Measured probe snapshot (Windows `.venv`, CPython 3.13, same machine/session):

- Baseline probe (pre-change snapshot, same loop structure):
  - `compute_sun`: `6.4689 µs/call`
  - `compute_moon`: `5.8760 µs/call`
  - `find_proper_suntime`: `6.2270 µs/call`
  - `find_proper_moontime`: `12.0292 µs/call`
  - `compute_visibilities_batch`: `23.4836 ms/call`
- Post-change probe (same probe):
  - `compute_sun`: `6.2658 µs/call` (~3.1% faster)
  - `compute_moon`: `5.8284 µs/call` (~0.8% faster)
  - `find_proper_suntime`: `6.2220 µs/call` (flat)
  - `find_proper_moontime`: `12.0750 µs/call` (within noise)
  - `compute_visibilities_batch`: `23.1195 ms/call` (~1.5% faster)
- New phase-4 probe (`examples/perf_phase4_probe.py --json`) additionally reports:
  - `find_utc_offset_varying_days_per_call_us`: `1291.26`
  - `find_utc_offset_warm_cache_per_call_us`: `0.617`

## Assumptions and Defaults

1. The initial audit pass was review/documentation-only; subsequent implementation updates are tracked in the phase execution status sections above.
2. Public API stability is preferred unless a correctness defect requires contract change.
3. Numeric correctness takes precedence over micro-optimizations.
4. Performance recommendations are prioritized by practical user impact (not purely theoretical gains).

## Completion Status

Status: **Complete for requested scope**

- Core Python and native C code paths were reviewed comprehensively.
- Findings are prioritized and implementation-ready.
- A staged roadmap is provided to execute changes with controlled regression risk.
