---
description: "Scaffold a C++ GTest or Python KratosUnittest test file for a feature in a Kratos application."
name: "Scaffold Test"
argument-hint: "Feature and application, e.g. 'ComputeVolumetricStrainProcess in StructuralMechanicsApplication (C++)'"
agent: agent
tools: [read, edit, search]
---

Generate a test file for a feature in this repository.

## Inputs

The user has provided: **{{input}}**

Parse:
- **FeatureName** — the class or process being tested (e.g. `ComputeVolumetricStrainProcess`)
- **ApplicationName** — target application folder under `applications/`
- **Language** — `C++` or `Python` (default: `C++` if not specified)

## Steps

1. **Read conventions** — load [testing instructions](../instructions/testing.instructions.md).

2. **Explore existing tests** in `applications/<ApplicationName>/tests/` to match the naming, fast-suite header, and suite registration style.

3. **For C++ tests**, create `applications/<ApplicationName>/tests/cpp_tests/test_<feature_name>.cpp`:
   - Include `testing/testing.h` and the application fast-suite header (e.g. `<app>_fast_suite.h`)
   - Use `KRATOS_TEST_CASE_IN_SUITE(...)` macro
   - Wrap in `namespace Kratos::Testing`
   - Cover at least: a happy-path case, a boundary/edge case, and a `Check()` validation case
   - Use `KRATOS_EXPECT_NEAR` (or appropriate `KRATOS_EXPECT_*` macro) for float comparisons

4. **For Python tests**, create `applications/<ApplicationName>/tests/test_<feature_name>.py`:
   - Import `KratosMultiphysics.KratosUnittest as KratosUnittest` and relevant Kratos modules
   - Class inheriting from `KratosUnittest.TestCase`
   - `setUp` creating a `Model` and `ModelPart`
   - At least one `test_` method per logical behaviour

5. **Summarise** the tests created and remind the user to add the new test file to the application suite runner (`test_<ApplicationName>.py`) if it is not auto-discovered.
