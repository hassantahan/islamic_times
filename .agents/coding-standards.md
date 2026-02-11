# Coding Standards

## General Principles
- Keep modifications minimal and task-focused.
- Match existing naming, module boundaries, and project idioms.
- Preserve behavior unless a task explicitly asks for a behavior change.

## Type Hints (Required)
- All new or modified Python functions and methods must include parameter and return type hints.
- New module-level constants should be explicitly typed when practical.
- Prefer specific container types over ambiguous untyped containers.
- Avoid `# type: ignore` where possible.
- If `# type: ignore` is necessary, keep scope narrow and add a short rationale comment.

## Docstrings (Required)
- Required for all public classes, methods, and functions.
- Required for non-trivial internal helpers.
- Include:
  - Purpose
  - Parameter meaning and units where relevant
  - Return value meaning and units where relevant
  - Exceptions raised for invalid input

## Comments (Required When Logic Is Non-Obvious)
- Add comments for algorithmic intent, domain assumptions, and edge-case handling.
- Do not add comments for obvious line-by-line operations.
- For astronomy/prayer calculations, explain why a branch/formula is used when not self-evident.

## Validation and Error Handling
- Keep input validation explicit and fail early.
- Use `TypeError` for wrong types and `ValueError` for invalid ranges/values.
- Keep error messages specific and actionable.

## Python Style Alignment
- Prefer readable explicit logic over compact clever constructs.
- Keep function responsibilities narrow and composable.
- Preserve existing public naming and data-shape conventions unless explicitly changing API.

## C-Extension Style Alignment
- Keep C and Python interfaces synchronized (names, argument order, expected units).
- When changing C calculations, update related comments and assumptions in code.
- Keep header/source contracts (`src/native/include/` and `src/native/`) in sync.
- Prefer deterministic behavior and explicit boundary checks.

## Domain and Terminology Consistency
- Preserve established Islamic astronomy/prayer terminology used by the project.
- Maintain unit consistency across APIs, dataclasses, and output formatting.
