---
name: scaffold-application
description: "Scaffold a complete new Kratos application. Use when creating a new application folder under applications/. Delegates to the repository's own application generator instead of hand-writing files."
argument-hint: "Application name in PascalCase, e.g. 'MyNewApplication'"
---

# Scaffold Application

## When to Use

Use this skill when the user asks to create a brand-new Kratos application under `applications/`.

## Procedure

1. **Do not hand-write the skeleton.** `kratos/python_scripts/application_generator/` already
   generates it: `createApplication.py` (and `laplacian_application_example.py`) show the
   `ApplicationGenerator` API — `AddVariables(...)`, `AddElements(...)`, `AddConditions(...)`,
   then `.Generate()`.

2. Copy one of those scripts, set the application name and the requested variables/elements/
   conditions, and run it (`python <your_script>.py`) from
   `kratos/python_scripts/application_generator/`. This produces the full
   `applications/<ApplicationName>/` tree: CMakeLists.txt, application header/source, variables
   header, pybind module, and test suite stub.

3. **Register the application for the build.** Applications are *not* wired into the root
   `CMakeLists.txt` — it only loops over the `KRATOS_APPLICATIONS` environment variable
   (`CMakeLists.txt:713-730`) and calls `add_subdirectory` for each entry. Add the new
   application's path to `KRATOS_APPLICATIONS` in `build/configure.sh` (or `.bat`).

4. Review the generated files against a small existing application (e.g.
   `ConvectionDiffusionApplication`) if anything looks off, and fill in real
   Elements/Conditions/Processes in place of the generator's placeholders.

5. **Summary** — report the files generated, the `KRATOS_APPLICATIONS` entry added, and remind
   the user to rebuild via `scripts/standard_configure.*` and to add the application to the CI
   JSON files under `.github/workflows/` if it needs CI coverage.
