---
type: C++ Function
title: BasicTestBuilderAndSolverDisplacement
resource: tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver.cpp#L73-L147
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void BasicTestBuilderAndSolverDisplacement( ModelPart& rModelPart, const bool WithConstraint = false, const bool AdditionalNode = false, const bool InvertRoleAdditionalNode = false )`

# Calls

- [CreateNewProperties](../../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [Fix](../../../../../../functions/includes/node/Fix.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/KRATOS_TEST_CASE_IN_SUITE.md)