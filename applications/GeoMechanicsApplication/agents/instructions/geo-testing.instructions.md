---
description: "Use when writing, reviewing or running tests for applications/GeoMechanicsApplication. Covers the GeoMechanics C++ test suites (with and without kernel), test helpers, tolerances, the Python suites in test_GeoMechanicsApplication.py, test data layout and JSON formatting."
applyTo: ["applications/GeoMechanicsApplication/tests/**", "applications/GeoMechanicsApplication/test_setup_utilities/**", "applications/GeoMechanicsApplication/benchmarks/**"]
---

# Testing conventions: GeoMechanicsApplication

The general rules in `/agents/instructions/testing.instructions.md` apply. This file adds what
is specific to this application.

## C++ unit tests

Location and build:

- Tests live in `tests/cpp_tests/<folder>/test_<name>.cpp`, where `<folder>` mirrors the source
  folder (`custom_elements`, `custom_utilities`, `custom_processes`, ...).
- `CMakeLists.txt` collects every `tests/*.cpp` with `GLOB_RECURSE` into the
  `KratosGeoMechanicsCoreTest` executable, so no CMake change is needed for a new test file.
- Every file starts with the GeoMechanics header block (see `geo-cpp.instructions.md`) and puts
  its tests in `namespace Kratos::Testing`.

Choose the suite (fixture):

| Suite | Header | Use when |
|-------|--------|----------|
| `KratosGeoMechanicsFastSuiteWithoutKernel` | `tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h` | Default. The test builds its own objects and only needs the variables registered in `geo_mechanics_fast_suite_without_kernel.cpp`. About 500 of the test cases use it |
| `KratosGeoMechanicsFastSuite` | `tests/cpp_tests/geo_mechanics_fast_suite.h` | The test needs the Kratos kernel with the GeoMechanics and LinearSolvers applications loaded (registered elements, conditions, laws, processes, Python-free model setup through the registry) |
| `KratosGeoMechanicsIntegrationSuite` | `tests/cpp_tests/geo_mechanics_fast_suite.h` | Slow tests that run a whole analysis. Excluded by `run_cpp_unit_tests.py --onlyFastSuite` |

If a test in the without-kernel suite needs a variable that is not registered there, add it to
`geo_mechanics_fast_suite_without_kernel.cpp` instead of moving the test to the kernel suite.

Macros and names:

- Use `KRATOS_TEST_CASE_IN_SUITE(TestName, Suite)` or `TEST_F(Suite, TestName)`; both are used.
  Use `TEST_P` with `INSTANTIATE_TEST_SUITE_P` for the same check over a table of inputs (see
  `custom_retention/test_van_genuchten_law.cpp`).
- Test names are `PascalCase` and say the subject and the expected behaviour, often as
  `ClassName_Behaviour`, for example `VanGenuchtenLawChecksInputParameters` or
  `MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureZone`.
- Prefer `KRATOS_EXPECT_*` from `includes/expect.h`. GTest macros (`EXPECT_NEAR`,
  `EXPECT_DOUBLE_EQ`) are accepted when no Kratos macro fits.
- Test error paths with `KRATOS_EXPECT_EXCEPTION_IS_THROWN(statement, "expected message")`
  and check the full message text.

Helpers (use them before writing new ones):

- `tests/cpp_tests/test_utilities.h`: `Defaults::absolute_tolerance` (1e-12) and
  `Defaults::relative_tolerance` (1e-6), `ExpectPointsAreNear`, and block checks for U-Pw
  systems (`AssertLHSMatrixBlocksAreNear`, `AssertUUBlockMatrixIsNear`,
  `AssertRHSVectorBlocksAreNear`, ...).
- `test_setup_utilities/model_setup_utilities.h` (`ModelSetupUtilities`) and
  `test_setup_utilities/element_setup_utilities.hpp` (`ElementSetupUtilities`) create model
  parts, nodes and elements. They are a separate library shared with `benchmarks/`.
- Stubs in `tests/cpp_tests/` (`stub_constitutive_law.h`, `stub_linear_elastic_law.h`,
  `stub_input_utility.h`, `stub_time_loop_executor.h`, `flow_stubs.h`, ...) replace heavy
  dependencies.

Expected values must come from an independent source (hand calculation, analytical solution or
a documented reference run), not from copying the current output of the code under test.

## Python tests

- Test modules are `tests/test_<topic>.py` (a few older ones have no `test_` prefix). The class
  name is `KratosGeoMechanics<Topic>Tests` and inherits `KratosUnittest.TestCase`, or
  `KratosGeoUnittest.TestCase` from `tests/KratosGeoUnittest.py` when it needs
  `assert_nodal_values_at_time` / `assert_integration_point_tensors` on GiD output.
- Each test method has a docstring with the physical set-up: geometry, material values with
  units, loading per stage, and the expected result.
- Use `tests/test_helper.py` to run cases and read results: `get_file_path`, `run_kratos`,
  `make_geomechanics_analysis`, `get_displacement`, `get_water_pressure`,
  `get_cauchy_stress_tensor`, `get_on_integration_points`, `calculate_delta`, ...
  Multi-stage cases use `KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages`.
  GiD results are read with `GiDOutputFileReader`.
- Input data of a case lives in `tests/<case_name>/` with `ProjectParameters*.json`,
  `MaterialParameters*.json`, the `.mdpa` mesh, and a `README.md` that describes the case.
  Shared files between variants go in a `common/` (or `Common/`) sub-folder.
- Register every new test class in `tests/test_GeoMechanicsApplication.py`:
  `small_test_cases` for short regression tests, `night_test_cases` for tests that take minutes,
  `valid_test_cases` for long validation runs. The nightly suite also runs all small tests.
- Plots: some tests write plots when the environment variable `KRATOS_GEO_MAKE_TEST_PLOTS` is
  `ON` (see `tests/README.md`). Plot code must not run otherwise.

## JSON input files

`exceptions_for_json_formatting.txt` lists the JSON files that are exempt from the JSON
formatting check. The checker itself is not in this repository. Write new JSON files with a
consistent 4-space indentation and do not add them to the exceptions list. When you reformat a
listed file, remove it from the list.

## Benchmarks

Google Benchmark files live in `benchmarks/` (for example `upw_diff_order_element_benchmark.cpp`).
Each `.cpp` becomes its own executable when `KRATOS_BUILD_BENCHMARK=ON`. Use the
`test_setup_utilities` library to build the model.

