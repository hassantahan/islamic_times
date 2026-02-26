# Coding Standards

## General Principles

- Keep modifications scoped to the request.
- Preserve public behavior unless a task explicitly calls for behavior change.
- Prefer readable, explicit code over compact but opaque logic.

## Python Requirements

- New and modified functions must include parameter and return type hints.
- Public and non-trivial functions require docstrings.
- Use NumPy-style sections in new/updated docstrings (`Parameters`, `Returns`, `Raises`, `Notes`) where applicable.
- Avoid broad `except Exception` blocks; translate only expected domain failures.
- Avoid `# type: ignore`; if required, keep scope narrow and include a short rationale.

## C Extension Requirements

- Keep headers and implementations synchronized (`src/native/include/*.h` and `src/native/*.c`).
- Keep Python wrapper contracts aligned with native contracts (argument order, units, error behavior).
- Add concise contract comments for non-obvious functions (inputs, units, status/error behavior).
- Prefer bounded loops and explicit failure paths over open-ended convergence loops.

## Comments and Documentation

- Add comments only for non-obvious intent, assumptions, or edge-case handling.
- Remove stale TODOs or replace them with clear issue references and scope.
- Keep public docs in `README.md`, `CONTRIBUTING.md`, and `docs/` user-focused.
- Do not add internal work-log documents to the repository unless explicitly requested.

## Terminology and Units

- Maintain established Islamic astronomy and prayer-time terminology.
- Keep units consistent across Python and C boundaries:
  - Angles: decimal degrees (unless explicitly documented otherwise)
  - UTC offsets: hours
  - Distances: kilometers
