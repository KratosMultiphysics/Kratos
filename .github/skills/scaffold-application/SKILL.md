---
name: scaffold-application
description: "Scaffold a complete new Kratos application. Use when creating a new application folder under applications/, including all required files: application header/source, CMakeLists.txt, variable registration, pybind11 module, and test suite stub."
argument-hint: "Application name in PascalCase, e.g. 'MyNewApplication'"
---

# Scaffold Application

## When to Use

Use this skill when the user asks to create a brand-new Kratos application under `applications/`.

## Procedure

### 1. Parse Inputs

Derive names from the user-supplied `ApplicationName` (PascalCase):

| Token | Example |
|-------|---------|
| `ApplicationName` | `MyNewApplication` |
| `application_name` | `my_new_application` |
| `APP_NAME` | `MY_NEW_APPLICATION` |
| `AppShort` (remove "Application") | `MyNew` |

### 2. Read Conventions

Load the following before generating any file:
- [cpp-conventions instructions](../../instructions/cpp-conventions.instructions.md)
- [cmake instructions](../../instructions/cmake.instructions.md)
- [python-conventions instructions](../../instructions/python-conventions.instructions.md)
- [testing instructions](../../instructions/testing.instructions.md)

### 3. Explore an Existing Application for Reference

Search `applications/` for a small existing application (e.g. `ConvectionDiffusionApplication`) and read its structure to match the exact style used in this repo.

### 4. Generate Files

Create the following files under `applications/<ApplicationName>/`:

```
applications/<ApplicationName>/
├── CMakeLists.txt
├── <ApplicationName>.py                    # Python __init__ for the application
├── <application_name>_application.h
├── <application_name>_application.cpp
├── <application_name>_application_variables.h
├── <application_name>_application_variables.cpp
├── custom_elements/          (empty, add .gitkeep)
├── custom_conditions/        (empty, add .gitkeep)
├── custom_processes/         (empty, add .gitkeep)
├── custom_utilities/         (empty, add .gitkeep)
├── custom_constitutive/      (empty, add .gitkeep)
├── python_scripts/
│   └── __init__.py
├── tests/
│   ├── __init__.py
│   └── test_<ApplicationName>.py   (minimal suite runner: small, nightly, validation, all)
└── custom_python/
    ├── <application_name>_python_application.cpp   (PYBIND11_MODULE entry)
    └── add_custom_processes_to_python.cpp           (empty binding stub)
```

**Key content rules:**
- `<application_name>_application.h` — declares `class <ApplicationName> : public KratosApplication`, includes `KRATOS_CLASS_POINTER_DEFINITION`
- `<application_name>_application.cpp` — implements `RegisterComponents()` and the factory entry point; calls `KRATOS_REGISTER_VARIABLE` for any new variables
- `<application_name>_application_variables.h` — declares new `KRATOS_DEFINE_APPLICATION_VARIABLE` entries
- `CMakeLists.txt` — follows the standard Kratos CMake pattern (see [cmake instructions](../../instructions/cmake.instructions.md))
- `custom_python/<application_name>_python_application.cpp` — `PYBIND11_MODULE` entry that calls all `AddCustom*ToPython` functions and registers variables
- `python_scripts/__init__.py` — imports the compiled pybind module

### 5. Update Root CMakeLists.txt

Check the root `CMakeLists.txt` for the pattern used to include other applications (typically an `option(KRATOS_BUILD_<APP> ...)` + `if(KRATOS_BUILD_<APP>) add_subdirectory(...) endif()` block) and add the new application following the same pattern.

### 6. Summary

Report:
- All files created
- The option/`add_subdirectory` lines added to root `CMakeLists.txt`
- Next steps: add Elements/Conditions/Processes, rebuild with `scripts/standard_configure.*`, and add the application to the CI JSON files in `.github/workflows/`
