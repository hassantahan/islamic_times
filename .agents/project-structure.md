# Project Structure

## Top-Level Layout
- `src/islamic_times/`: Main Python package (src-layout).
- `src/native/`: Native C extension implementation files.
- `src/native/include/`: C headers used by the native extension sources.
- `tests/`: Automated test suite for package behavior.
- `examples/`: Runnable usage examples.
- `mapper.py`: New moon visibility mapping and plotting workflow.
- `setup.py`: Extension build configuration.
- `pyproject.toml`: Build backend, project metadata, extras, and tool configuration.

## Repository Tree (With Short Descriptions)
```text
islamic_times/
├── .agents/                              # Agent-specific project guidance
│   ├── README.md                         # Index and read order for agent docs
│   ├── coding-standards.md               # Required style, typing, docs, comments
│   ├── project-structure.md              # Repo structure and ownership map
│   └── workflow-checklist.md             # Implementation and validation checklist
├── src/                                  # Source root (Python + native)
│   ├── islamic_times/                    # Python package source (src-layout)
│   │   ├── __init__.py                   # Package init and version fallback
│   │   ├── calculation_equations.py      # Pure-Python mathematical utilities
│   │   ├── islamic_times.py              # `ITLocation` public orchestration class
│   │   ├── it_dataclasses.py             # Typed dataclasses/value objects
│   │   ├── moon_equations.py             # Moon position/phase algorithms
│   │   ├── prayer_times.py               # Prayer-time computation logic
│   │   ├── sun_equations.py              # Sun position/event algorithms
│   │   └── time_equations.py             # Calendar/timezone/time conversions
│   └── native/                           # Native C extension implementation
│       ├── astro_core.c                  # Python C API module bindings
│       ├── c_calculation_equations.c     # C math helper implementations
│       ├── c_datetime.c                  # C datetime conversion implementations
│       ├── c_moon_equations.c            # C lunar calculation implementations
│       ├── c_sun_equations.c             # C solar calculation implementations
│       ├── c_time_equations.c            # C astronomical time implementations
│       ├── c_visibilities.c              # C crescent visibility implementations
│       └── include/                      # Headers for native extension modules
│           ├── astro_core.h              # Top-level C extension interface
│           ├── c_calculation_equations.h # Math helper declarations
│           ├── c_datetime.h              # Datetime conversion declarations
│           ├── c_moon_equations.h        # Lunar calculation declarations
│           ├── c_sun_equations.h         # Solar calculation declarations
│           ├── c_time_equations.h        # Astronomical time declarations
│           └── c_visibilities.h          # Crescent visibility declarations
├── tests/                                # Automated tests for library behavior
├── examples/demo.py                      # Small runnable API example
├── map_shp_files/                        # Shapefiles used by mapping pipeline
├── test_maps/                            # Generated/fixture visibility map outputs
├── mapper.py                             # CLI/script for visibility map generation
├── CONTRIBUTING.md                       # Contributor setup and validation guide
├── LICENSE                               # Project license
├── MANIFEST.in                           # Source distribution include rules
├── pyproject.toml                        # Build-system and tool configuration
├── README.md                             # Public package docs and usage examples
├── RELEASE_CHECKLIST.md                  # Release gating checklist
└── setup.py                              # Package + C-extension build setup
```

## Core Runtime Flow
1. `src/islamic_times/islamic_times.py` provides `ITLocation`, the primary orchestration class.
2. `ITLocation` validates observer/date inputs and coordinates computations.
3. Heavy astronomy math is delegated to `islamic_times.astro_core` (C extension).
4. Python modules (`prayer_times.py`, `moon_equations.py`, `sun_equations.py`, `time_equations.py`) provide domain logic and integration behavior.
5. Dataclasses and typed value objects live in `src/islamic_times/it_dataclasses.py`.

## Where To Make Changes
- Prayer-time algorithm changes: start in `src/islamic_times/prayer_times.py`; keep contracts compatible with `ITLocation`.
- Astronomical core logic changes:
  - Python layer behavior: `src/islamic_times/sun_equations.py`, `src/islamic_times/moon_equations.py`, `src/islamic_times/time_equations.py`, `src/islamic_times/calculation_equations.py`.
  - Performance-sensitive numeric core: `src/native/*.c` and matching `src/native/include/*.h`.
- Output shape/type changes: verify impacts in `src/islamic_times/it_dataclasses.py` and all callers.
- Mapping changes: `mapper.py` plus any data-flow assumptions for visibility outputs.

## API and Behavior Contracts To Preserve
- `ITLocation` is the public entrypoint and should remain stable unless explicitly requested.
- Units and semantics must remain clear and consistent:
  - Angles in decimal degrees unless otherwise specified.
  - UTC offsets in hours.
  - Datetime values with explicit timezone handling.
- Prayer outputs must preserve established ordering and meaning (Fajr, Sunrise, Zuhr, Asr, Sunset, Maghrib, Isha, Midnight).
- If changing behavior intentionally, document the compatibility impact in the change summary.
