---
name: geo-build-and-test
description: "Use when building GeoMechanicsApplication or running its C++ or Python tests. Extends the root build skill with the GeoMechanics dependencies, test executables, gtest filters and Python suites."
---

# Build and test GeoMechanicsApplication

Read `/agents/skills/build/SKILL.md` first. Build only through the personal copy of the
configure script (`build/configure.bat` or `build/configure.sh`), never with direct `cmake` calls.

## Configure

GeoMechanics needs two other applications. In `build/configure.bat`:

```bat
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\GeoMechanicsApplication;
```

In `build/configure.sh`:

```sh
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/GeoMechanicsApplication
```

Pass `-DKRATOS_BUILD_TESTING=ON` to build the C++ tests, and `-DKRATOS_BUILD_BENCHMARK=ON` for
`benchmarks/`. Do not enable MPI for this application; it is not supported.

Build targets produced by `applications/GeoMechanicsApplication/CMakeLists.txt`:

| Target | Contents |
|--------|----------|
| `KratosGeoMechanicsCore` | Application library |
| `KratosGeoMechanicsApplication` | Python module |
| `KratosGeoMechanicsTestSetupUtilities` | Shared test and benchmark helpers |
| `KratosGeoMechanicsCoreTest` | All C++ tests under `tests/` |

When CMake was already configured, a rebuild of one target is enough after a source change.
A new `.cpp` file needs a new configure run because sources are collected with `GLOB_RECURSE`.

## Run C++ tests

The executable is installed as `bin/<BuildType>/test/KratosGeoMechanicsCoreTest` (`.exe` on
Windows). The `PATH` (Windows) or `LD_LIBRARY_PATH` (Linux) must contain `bin/<BuildType>/libs`.

Run the narrowest filter first:

```powershell
$env:PATH = "C:\Checkouts\KratosProjects\yek\bin\Release\libs;$env:PATH"
& C:\Checkouts\KratosProjects\yek\bin\Release\test\KratosGeoMechanicsCoreTest.exe --gtest_filter="*VanGenuchten*"
```

Then the whole fast set (both fast suites, no integration tests):

```powershell
& ...\KratosGeoMechanicsCoreTest.exe --gtest_filter="KratosGeoMechanicsFastSuite*"
```

`tests/run_cpp_unit_tests.py --onlyFastSuite` does the same when started from `bin/<BuildType>`.

## Run Python tests

The Python test commands need this environment; without it, `import KratosMultiphysics` fails.
If you use your own shell aliases or wrapper scripts that set it up, those work too.

```powershell
$env:PYTHONPATH = "<kratos_root>\bin\Release"
$env:PATH = "<kratos_root>\bin\Release\libs;$env:PATH"
Set-Location <kratos_root>\applications\GeoMechanicsApplication\tests
```

On Linux, set `PYTHONPATH` and `LD_LIBRARY_PATH` in the same way. Run from
`applications/GeoMechanicsApplication/tests/` so the test modules and data are found.

- One test class or method:
  `python -m unittest test_dirichlet_release.KratosGeoMechanicsDirichletReleaseTests`
- One suite: `python test_GeoMechanicsApplication.py -l small` (also `nightly`, `validation`,
  `all`).

Set `KRATOS_GEO_MAKE_TEST_PLOTS=ON` only when you want the result plots that some tests make.

## Before you report success

1. The changed C++ files are formatted with `applications/GeoMechanicsApplication/.clang-format`.
2. The new or changed C++ tests pass with a filter that selects them.
3. The full fast C++ set passes.
4. The Python tests that use the changed code pass; for a change in a shared element, law or
   process, the `small` suite passes.
