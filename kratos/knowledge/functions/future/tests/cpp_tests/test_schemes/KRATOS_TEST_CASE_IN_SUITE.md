---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: future/tests/cpp_tests/test_schemes.cpp#L37-L89
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/SolvingStrategiesTestUtilities/SetUpTestModelPart1D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/SolvingStrategiesTestUtilities/SetUpTestModelPart2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildWithSafeAssemble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/BuildWithThreadLocal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/BuildWithLocalAllocation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(StaticSchemeBuild1D, KratosCoreFastSuite)`

# Calls

- [SetUpTestModelPart1D](../../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/SolvingStrategiesTestUtilities/SetUpTestModelPart1D.md)
- [pGetVector](../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)
- [SetUpTestModelPart2D](../../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/SolvingStrategiesTestUtilities/SetUpTestModelPart2D.md)
- [BuildWithSafeAssemble](../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildWithSafeAssemble.md)
- [BuildWithThreadLocal](../../../../../functions/future/solving_strategies/schemes/implicit_scheme/BuildWithThreadLocal.md)
- [BuildWithLocalAllocation](../../../../../functions/future/solving_strategies/schemes/implicit_scheme/BuildWithLocalAllocation.md)