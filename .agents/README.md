# .agents Authoring Guide

## What This Folder Is For
This folder defines project-specific instructions for coding agents working on `islamic_times`.

## Project Snapshot
`islamic_times` is a Python package with a C extension for Islamic astronomical calculations, prayer times, Hijri conversion, and moon visibility analysis.

## Read Order
1. `project-structure.md`
2. `coding-standards.md`
3. `workflow-checklist.md`

## Scope
These guides apply to:
- Python package code in `src/islamic_times/`
- C extension sources in `src/native/`
- C headers in `src/native/include/`
- Supporting scripts such as `mapper.py` and `examples/demo.py`

## Priority Rules
1. Follow existing project patterns unless a task explicitly asks for a refactor.
2. Preserve public behavior and data semantics unless the task explicitly changes them.
3. Keep changes targeted, typed, and documented.
