---
name: ci-cd
description: "Use when modifying or reviewing files under .github/workflows/ in Kratos Multiphysics. Not needed for everyday code changes."
---

# CI/CD Conventions — Kratos Multiphysics

Only relevant when a task explicitly touches `.github/workflows/`. Otherwise leave CI alone.

## Platform

This project uses **GitHub Actions** exclusively — do **not** suggest GitLab CI, Travis CI, or
any other CI system. Secrets/env vars are managed via GitHub Actions secrets, never hardcoded.

## Key Workflows

| File | Trigger | Purpose |
|------|---------|---------|
| `ci.yml` | PRs to `master`, `workflow_dispatch` | Main PR check — builds and tests on Ubuntu matrix |
| `nightly_build.yml` | Scheduled nightly | Broader coverage including Windows and Rocky Linux |

`ci.yml` first computes changed files/applications (to skip unaffected builds), then runs an
Ubuntu build/test matrix; the nightly workflow adds Windows and Rocky Linux jobs.

## Application Selection Files

Which applications get built per CI environment is controlled by JSON files in
`.github/workflows/`: `ci_apps_linux.json`, `ci_apps_windows.json`, `ci_apps_rocky.json`,
`ci_apps_intel.json`. Add a new application to the relevant file(s) when it needs CI coverage.

## Rules

- Preserve existing job naming and `on:` / `strategy:` / `concurrency:` style.
- Keep the `concurrency` group (`ci-${{ github.head_ref }}`) and `cancel-in-progress: true`.
- Only change CI behavior when explicitly requested; do not refactor unrelated jobs.
