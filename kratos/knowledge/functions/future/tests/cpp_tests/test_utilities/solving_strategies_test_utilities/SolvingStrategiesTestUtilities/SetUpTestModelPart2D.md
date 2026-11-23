---
type: C++ Method
title: SetUpTestModelPart2D
resource: future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities.h#L320-L354
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/tests/cpp_tests/test_schemes/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void SetUpTestModelPart2D( const std::size_t NumElemsX, const std::size_t NumElemsY, const double ElemSizeX, const double ElemSizeY, ModelPart& rModelPart)`

# Calls

- [X0](../../../../../../../functions/includes/node/X0.md)
- [Y0](../../../../../../../functions/includes/node/Y0.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../../../functions/future/tests/cpp_tests/test_schemes/KRATOS_TEST_CASE_IN_SUITE.md)