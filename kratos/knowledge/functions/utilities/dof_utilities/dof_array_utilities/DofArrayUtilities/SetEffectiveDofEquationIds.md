---
type: C++ Method
title: SetEffectiveDofEquationIds
resource: utilities/dof_utilities/dof_array_utilities.cpp#L181-L208
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/SetEffectiveEquationId
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

`void DofArrayUtilities::SetEffectiveDofEquationIds( const DofsArrayType& rDofArray, DofsArrayType& rEffectiveDofArray)`

# Calls

- [SetEffectiveEquationId](../../../../../functions/includes/dof/Dof/SetEffectiveEquationId.md)

# Called by

- [SetUpSystemIds](../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpSystemIds.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE.md)