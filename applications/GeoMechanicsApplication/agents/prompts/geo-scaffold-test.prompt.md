---
description: "Create a C++ unit test or a Python integration test for a feature in GeoMechanicsApplication, following the GeoMechanics test conventions."
name: "Scaffold GeoMechanics Test"
argument-hint: "Feature and language, e.g. 'VanGenuchtenLaw (C++)' or 'ApplyK0ProcedureProcess on a layered column (Python)'"
agent: agent
tools: [read, edit, search]
---

Create a test for a feature in `applications/GeoMechanicsApplication/`.

## Input

The user has provided: **{{input}}**

From it, find:
- **Feature**: the class, process or behaviour under test
- **Language**: `C++` or `Python` (default `C++`)

## Steps

1. Read `applications/GeoMechanicsApplication/agents/instructions/geo-testing.instructions.md`
   and the source of the feature.
2. Find existing tests of the same folder in `tests/cpp_tests/<folder>/` (C++) or `tests/`
   (Python) and match their style. If a test file for the feature already exists, add to it.
3. **C++**: create `tests/cpp_tests/<mirrored folder>/test_<feature_name>.cpp` with:
   - the GeoMechanics header block (ask the user for the author name if it is unknown);
   - `#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"` (or
     `geo_mechanics_fast_suite.h` when the kernel is needed) and `includes/expect.h`;
   - `namespace Kratos::Testing`;
   - tests named `<Feature>_<Behaviour>` for the main result, one edge case, and each `Check()`
     error (`KRATOS_EXPECT_EXCEPTION_IS_THROWN` with the full message);
   - helpers from `tests/cpp_tests/test_utilities.h` and `test_setup_utilities/`, tolerances
     from `Defaults::`;
   - expected values from a hand calculation or analytical solution, with a short comment that
     gives the formula or source.
4. **Python**: create the case folder `tests/<case_name>/` with `ProjectParameters.json`,
   `MaterialParameters.json` (4-space indent), the mesh and a `README.md` that describes the
   set-up. Create or extend `tests/test_<topic>.py` with a class
   `KratosGeoMechanics<Topic>Tests(KratosUnittest.TestCase)`, a docstring per test with the
   physical set-up and units, and `test_helper` functions to run and read results. Add the
   class to the right list in `tests/test_GeoMechanicsApplication.py`.
5. Report the files created and the command that runs only the new tests (see
   `agents/skills/geo-build-and-test/SKILL.md`).
