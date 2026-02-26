# Project Structure

## Top-Level Layout

- `src/islamic_times/`: Main Python package source.
- `src/islamic_times/mapper/`: Mapper pipeline package (CLI, compute, rendering, geodata, profiling).
- `src/islamic_times/_legacy_py_impl/`: Deprecated pure-Python algorithm implementations kept behind compatibility wrappers.
- `src/native/`: Native C extension implementation.
- `src/native/include/`: Native C headers.
- `tests/`: Automated test suite.
- `docs/`: Public repository documentation only.
- `.agents/`: Internal agent instructions.

## Key Runtime Ownership

1. `src/islamic_times/islamic_times.py`
   Public `ITLocation` orchestration API.
2. `src/islamic_times/prayer_times.py`
   Prayer-time logic, including extreme-latitude handling.
3. `src/islamic_times/sun_equations.py`, `src/islamic_times/moon_equations.py`, `src/islamic_times/time_equations.py`
   Domain integration around native astronomy kernels.
4. `src/native/astro_core.c`
   Python C-extension entrypoint and wrapper surface.
5. `src/native/c_event_solver.c`
   Shared event-search and interpolation helpers used by sun/moon paths.
6. `src/islamic_times/mapper/`
   Current mapper workflow implementation.
7. `mapper.py`
   Legacy compatibility shim delegating to `islamic_times.mapper`.

## Edit Targets by Change Type

- Prayer behavior changes:
  `src/islamic_times/prayer_times.py` and `src/islamic_times/islamic_times.py`
- Sun/moon/time algorithm or wrapper changes:
  `src/native/*.c`, `src/native/include/*.h`, plus matching Python integration modules
- Mapper behavior/performance changes:
  `src/islamic_times/mapper/*.py` and mapper-focused tests
- API surface/value-shape changes:
  `src/islamic_times/it_dataclasses.py`, `src/islamic_times/islamic_times.py`, tests, and public docs

## Contract Guardrails

- Keep `ITLocation` behavior stable unless explicitly changing public contract.
- Keep criterion and method enums consistent across Python and native layers.
- If legacy compatibility behavior changes, document migration impact in public docs.
