---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/strategies/convergence_criteria/test_residual_criteria.cpp#L51-L113
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/strategies/convergence_criteria/test_residual_criteria/GenerateTestResidualCriteriaModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/SetEquationId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(ResidualCriteria, KratosCoreFastSuite)`

# Calls

- [GenerateTestResidualCriteriaModelPart](../../../../../../functions/tests/cpp_tests/strategies/convergence_criteria/test_residual_criteria/GenerateTestResidualCriteriaModelPart.md)
- [SetEquationId](../../../../../../functions/includes/dof/Dof/SetEquationId.md)
- [Fix](../../../../../../functions/includes/node/Fix.md)