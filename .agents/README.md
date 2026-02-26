# .agents Authoring Guide

## Purpose

This folder contains internal guidance for coding agents working on `islamic_times`.
It is intentionally repository-tracked, but it is not part of the public package API docs.

## Read Order

1. `project-structure.md`
2. `coding-standards.md`
3. `workflow-checklist.md`

## Scope

These guides apply to:

- Python package code in `src/islamic_times/`
- Native extension code in `src/native/` and `src/native/include/`
- Tests in `tests/`
- Mapping pipeline code in `src/islamic_times/mapper/` and compatibility shim `mapper.py`

## Documentation Boundary

Public repository docs should stay in `README.md`, `CONTRIBUTING.md`, and `docs/`.
Internal agent guidance stays in `.agents/`.
