---
name: geo-add-process
description: "Use when adding a new C++ process to GeoMechanicsApplication. Covers the header and source, Python binding, Python Factory wrapper, optional registration in the C++ settlement workflow, unit tests, an integration test and the README section."
argument-hint: "Process name in PascalCase, e.g. 'ApplyInitialPorosityProcess'"
---

# Add a process to GeoMechanicsApplication

Read `AGENTS.md`, `geo-cpp.instructions.md`, `geo-python.instructions.md`,
`geo-testing.instructions.md` and `custom_processes/README.md` (all in
`applications/GeoMechanicsApplication/`) before you start.

Paths below are relative to `applications/GeoMechanicsApplication/`. `<Name>` is the
`PascalCase` process name, `<name>` its `snake_case` form.

## 1. Choose a reference

Use `custom_processes/apply_k0_procedure_process.h/.cpp` as the main reference. For a process
that only sets nodal values, also look at `geo_apply_constant_scalar_value_process.h/.cpp`.
If a core Kratos process with the same name exists, prefix the new name with `Geo`.

## 2. Header and source

Create `custom_processes/<name>.h` and `custom_processes/<name>.cpp`:

- GeoMechanics header block with the real author names, `#pragma once`.
- `class KRATOS_API(GEO_MECHANICS_APPLICATION) <Name> : public Process` with
  `KRATOS_CLASS_POINTER_DEFINITION(<Name>);`.
- Inspect the process requirements and nearby implementations to choose the constructor
  signature and model/model-part access pattern. When the process uses `model_part_name` or
  `model_part_name_list`, prefer a `<Name>(Model& rModel, Parameters Settings)` constructor
  that gets the model parts with
  `ProcessUtilities::GetModelPartsFromSettings(rModel, mSettings, Info())` from
  `custom_utilities/process_utilities.h`.
- Override only the stage methods the process needs, plus `int Check() override` when it has
  input to validate and `[[nodiscard]] std::string Info() const override`.
- Validate input with `KRATOS_ERROR_IF` / `CheckProperties` and messages that name the
  parameter and the model part.
- Format both files with `.clang-format`.

No CMake change is needed; `custom_processes/*.cpp` is globbed.

## 3. Python binding and wrapper

- In `custom_python/add_custom_processes_to_python.cpp`, include the header and add, next to
  the other processes:

  ```cpp
  py::class_<<Name>, <Name>::Pointer, Process>(m, "<Name>", "<one-line description>")
      .def(py::init<Model&, Parameters>());
  ```

- Create `python_scripts/<name>.py` with only a `Factory(settings, model)` function, as in
  `python_scripts/apply_k0_procedure_process.py`.

## 4. C++ settlement workflow (only if needed)

If the process must run in the C++ route (`custom_workflows/dgeosettlement.cpp`), register it
there with `mProcessFactory->AddCreator("<Name>", MakeCreatorWithModelFor<<Name>>());` and add or
extend a `*CppRoute` Python test.

## 5. Tests

- C++: `tests/cpp_tests/custom_processes/test_<name>.cpp` (see
  `test_apply_k0_procedure_process.cpp`). Cover the process behaviour and relevant edge cases
  or invalid inputs (test `Check()` errors with `KRATOS_EXPECT_EXCEPTION_IS_THROWN`); do not
  duplicate coverage already provided by existing tests. Use
  `KratosGeoMechanicsFastSuiteWithoutKernel` unless the kernel is needed.
- Python: reuse existing case data where suitable. Otherwise add a case folder
  `tests/<case_name>/` with `ProjectParameters.json` (4-space indent),
  `MaterialParameters.json`, the mesh and a `README.md`. Add a test class in
  `tests/test_<topic>.py` and an entry in the right list of `tests/test_GeoMechanicsApplication.py`.

## 6. Documentation

Add a section to `custom_processes/README.md` and a link in its list of documented processes:
purpose, JSON parameters with units and defaults, and a JSON example.

## 7. Verify and report

Follow `agents/skills/geo-build-and-test/SKILL.md`: rebuild (new `.cpp` files need a configure
run), run the new C++ tests by filter, the fast C++ set, and the new Python test. Report the
files created, the registrations added, and the test results.

