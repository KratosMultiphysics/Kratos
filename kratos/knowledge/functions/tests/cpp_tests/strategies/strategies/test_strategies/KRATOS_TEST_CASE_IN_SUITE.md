---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/strategies/strategies/test_strategies.cpp#L146-L180
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/strategies/strategies/test_strategies/BasicTestStrategyDisplacement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/SetUpNonconvergedSolutionsFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/GetNonconvergedSolutions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(DisplacementLinearStrategy, KratosCoreFastSuite)`

# Calls

- [BasicTestStrategyDisplacement](../../../../../../functions/tests/cpp_tests/strategies/strategies/test_strategies/BasicTestStrategyDisplacement.md)
- [SetUpNonconvergedSolutionsFlag](../../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/SetUpNonconvergedSolutionsFlag.md)
- [GetNonconvergedSolutions](../../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/GetNonconvergedSolutions.md)