---
type: C++ Method
title: SetDofEquationIds
resource: utilities/dof_utilities/dof_array_utilities.cpp#L172-L179
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/SetEquationId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpSystemIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DofArrayUtilities::SetDofEquationIds(const DofsArrayType &rDofArray)`

# Calls

- [SetEquationId](../../../../../functions/includes/dof/Dof/SetEquationId.md)

# Called by

- [SetUpSystemIds](../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpSystemIds.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE.md)