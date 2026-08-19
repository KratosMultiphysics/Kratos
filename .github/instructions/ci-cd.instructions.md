---
description: "Use when modifying or reviewing the GitHub Actions CI/CD workflows in Kratos Multiphysics. Covers workflow structure, job naming, matrix strategy, and constraints."
applyTo: ".github/workflows/*.yml"
---

# CI/CD Conventions — Kratos Multiphysics

## Platform

This project uses **GitHub Actions** exclusively. Do **not** suggest GitLab CI, Travis CI, or any other CI system.

- Workflow definitions: `.github/workflows/` at the repository root.
- Secrets and environment variables are managed through GitHub Actions repository/organization secrets — never hardcode them in workflow files.

## Key Workflows

| File | Trigger | Purpose |
|------|---------|---------|
| `ci.yml` | PRs to `master`, `workflow_dispatch` | Main PR check — builds and tests on Ubuntu matrix |
| `nightly_build.yml` | Scheduled nightly | Broader coverage including Windows and Rocky Linux |

## `ci.yml` Structure

- **`changed-files` job**: determines which applications have changed (via `get_files_changed_in_pr.py`) to skip unaffected builds.
- **`ubuntu` job**: matrix of `build-type` (Custom, FullDebug) × `compiler` (gcc, clang); runs on `kratosmultiphysics/kratos-image-ci-ubuntu-22-04:latest`.
- **`windows` job** (nightly): MSVC builds on Windows runners.
- **`rocky` job** (nightly): gcc build on Rocky Linux container.

## Environment Variables in CI

| Variable | Purpose |
|----------|---------|
| `KRATOS_BUILD_TYPE` | Taken from matrix |
| `KRATOS_CI_CHANGED_FILES` | JSON list of changed files, used to skip unaffected apps |
| `KRATOS_CI_APPLICATIONS` | JSON file listing which apps to build (e.g. `ci_apps_linux.json`) |
| `KRATOS_CI_CORES` | Parallel compile jobs (typically 4 for hosted runners) |
| `OMPI_MCA_rmaps_base_oversubscribe` | Allow OpenMPI oversubscription in containers |

## Application Selection Files

The following JSON files in `.github/workflows/` define which applications are built per CI environment:

- `ci_apps_linux.json` — Ubuntu PR builds
- `ci_apps_windows.json` — Windows nightly builds
- `ci_apps_rocky.json` — Rocky Linux nightly builds
- `ci_apps_intel.json` — Intel LLVM builds

When adding a new application to CI, add it to the appropriate JSON file(s).

## Rules

- Preserve job naming and `on:` / `strategy:` / `concurrency:` style consistent with existing workflows.
- Keep the `concurrency` group (`ci-${{ github.head_ref }}`) and `cancel-in-progress: true` to avoid wasted runner time.
- Only change CI behavior when explicitly requested; do not refactor unrelated jobs.
- Do **not** suggest non-GitHub CI systems.
