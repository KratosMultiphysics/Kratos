---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/sources/test_dof.cpp#L26-L37
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FixDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FreeDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/HasReaction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/SetEquationId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/SetEffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/EffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(DofConstructorWtihoutVariableInVariablesList, KratosCoreFastSuite)`

# Calls

- [IsFree](../../../../../functions/includes/dof/Dof/IsFree.md)
- [FixDof](../../../../../functions/includes/dof/Dof/FixDof.md)
- [FreeDof](../../../../../functions/includes/dof/Dof/FreeDof.md)
- [GetVariable](../../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [HasReaction](../../../../../functions/includes/dof/Dof/HasReaction.md)
- [SetEquationId](../../../../../functions/includes/dof/Dof/SetEquationId.md)
- [SetEffectiveEquationId](../../../../../functions/includes/dof/Dof/SetEffectiveEquationId.md)
- [EffectiveEquationId](../../../../../functions/includes/dof/Dof/EffectiveEquationId.md)