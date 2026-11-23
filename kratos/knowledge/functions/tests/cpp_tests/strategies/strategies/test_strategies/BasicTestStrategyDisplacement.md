---
type: C++ Function
title: BasicTestStrategyDisplacement
resource: tests/cpp_tests/strategies/strategies/test_strategies.cpp#L94-L140
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/strategies/strategies/test_strategies/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline DofsArrayType BasicTestStrategyDisplacement( ModelPart& ModelPart, const ResidualType ThisResidualType )`

# Calls

- [CreateNewProperties](../../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../../functions/tests/cpp_tests/strategies/strategies/test_strategies/KRATOS_TEST_CASE_IN_SUITE.md)