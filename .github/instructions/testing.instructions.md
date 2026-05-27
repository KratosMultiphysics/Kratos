---
description: "Use when writing, reviewing, or running tests. Covers C++ (GTest/Kratos), Python (KratosUnittest), benchmarks (Google Benchmark), and MPI test conventions for Kratos Multiphysics."
applyTo: ["**/tests/**", "**/test_*.cpp", "**/test_*.py", "**/benchmark_*.cpp"]
---

# Testing Conventions — Kratos Multiphysics

## C++ Tests

- Framework: **Kratos wrappers over GTest** (`testing/testing.h`).
- Include the application fast-suite header for fixture setup.
- Preferred test macro: `KRATOS_TEST_CASE_IN_SUITE(TestName, SuiteName)`.
- Enclose tests in `namespace Kratos::Testing`.
- Assertion macros live in `includes/expect.h` — **prefer `KRATOS_EXPECT_*` over raw GTest macros**.
- Test files are named `test_<feature_name>.cpp` and live in `tests/cpp_tests/` within each application.
- Each application defines a fast-suite fixture (e.g., `KratosStructuralMechanicsFastSuite`) inheriting from `KratosCoreFastSuite`.

```cpp
#include "testing/testing.h"
#include "structural_mechanics_fast_suite.h"   // adjust per application
#include "containers/model.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(MyFeature_ShouldDoSomething, KratosStructuralMechanicsFastSuite)
{
    // Arrange
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("TestPart", 1);

    // Act
    // ...

    // Assert — prefer KRATOS_EXPECT_* macros
    KRATOS_EXPECT_NEAR(result, expected, 1e-10);
    KRATOS_EXPECT_EQ(count, expected_count);
    KRATOS_EXPECT_VECTOR_NEAR(result_vector, expected_vector, 1e-10);
    KRATOS_EXPECT_MATRIX_NEAR(result_matrix, expected_matrix, 1e-10);
}

} // namespace Kratos::Testing
```

| Macro | Purpose |
|-------|---------|
| `KRATOS_EXPECT_NEAR(a, b, tol)` | Float equality with tolerance |
| `KRATOS_EXPECT_EQ(a, b)` | Exact equality |
| `KRATOS_EXPECT_TRUE(cond)` | Boolean assertion |
| `KRATOS_EXPECT_VECTOR_NEAR(v1, v2, tol)` | Vector near equality |
| `KRATOS_EXPECT_MATRIX_NEAR(m1, m2, tol)` | Matrix near equality |
| `KRATOS_EXPECT_EXCEPTION_IS_THROWN(stmt, msg)` | Exception assertion |

Use `KRATOS_CHECK_*` macros (from `includes/checks.h`) in **production code** preconditions, not in tests.

## Python Tests

- Framework: **`KratosMultiphysics.KratosUnittest`** (wrapper over `unittest`).
- Test files are named `test_<feature_name>.py` and placed in `tests/` within the application.
- Test classes must inherit from `KratosUnittest.TestCase`.
- Every test method must start with `test_`.
- Use Kratos-specific helpers on `TestCase`:
  - `skipTestIfApplicationsNotAvailable(*apps)`
  - `assertVectorAlmostEqual(v1, v2, places=7)`
  - `assertMatrixAlmostEqual(m1, m2, places=7)`
- Application suite files (`test_<ApplicationName>.py`) assemble `small`, `nightly`, `validation`, and `all` suites.

```python
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestMyFeature(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("TestPart")

    def test_something(self):
        # Arrange
        # Act
        # Assert
        self.assertAlmostEqual(expected, actual, places=6)
        self.assertVectorAlmostEqual(expected_vec, actual_vec)

if __name__ == "__main__":
    KratosUnittest.main()
```

## C++ Benchmarks

- Framework: **Google Benchmark**.
- Benchmark files are named `benchmark_<feature_name>.cpp` and live in `benchmarks/` within the application or core.
- Always use `benchmark::State` and the `BENCHMARK()` macro.
- End the file with `BENCHMARK_MAIN()`.
- Built only when `KRATOS_BUILD_BENCHMARK=ON`.

```cpp
#include <benchmark/benchmark.h>
#include "my_feature.h"

static void BM_MyFeature(benchmark::State& state) {
    // setup outside the loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(/* ... */);
    }
}

BENCHMARK(BM_MyFeature)->Range(8, 8 << 10);
BENCHMARK_MAIN();
```

## MPI Tests

- Use the dedicated VS Code tasks `Run C++ Tests` / `Run Tests` with MPI variants.
- Keep `OMP_NUM_THREADS=1` consistent with existing tasks and CI configuration.
- MPI test runner scripts are present under `bin/<BuildType>/KratosMultiphysics/testing/`.

## Running Tests Locally

Prefer the VS Code tasks rather than ad-hoc commands:

| Task | Purpose |
|------|---------|
| `Run Tests` | All Python test suites |
| `Run C++ Tests` | All C++ GTest suites |
| `Run C++ Test Suite` | Specific GTest executable |
| `Run C++ Test Suite Filtered` | Filtered GTest run |
| `Run Current Benchmark file to JSON` | Single benchmark executable |
