---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_dof_array_utilities.cpp#L58-L72
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/utilities/test_dof_array_utilities/SetTestDofArrayModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpEffectiveDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetDofEquationIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetEffectiveDofEquationIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/EffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(SetUpDofArray, KratosCoreFastSuite)`

# Calls

- [SetTestDofArrayModelPart](../../../../../functions/tests/cpp_tests/utilities/test_dof_array_utilities/SetTestDofArrayModelPart.md)
- [GetVariable](../../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [SetUpEffectiveDofArray](../../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpEffectiveDofArray.md)
- [SetDofEquationIds](../../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetDofEquationIds.md)
- [SetEffectiveDofEquationIds](../../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetEffectiveDofEquationIds.md)
- [EffectiveEquationId](../../../../../functions/includes/dof/Dof/EffectiveEquationId.md)