---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/sources/test_node.cpp#L28-L50
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/HasReaction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/GetDofPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(NodeAssignOperator, KratosCoreFastSuite)`

# Calls

- [Fix](../../../../../functions/includes/node/Fix.md)
- [IsFree](../../../../../functions/includes/dof/Dof/IsFree.md)
- [HasReaction](../../../../../functions/includes/dof/Dof/HasReaction.md)
- [GetDofPosition](../../../../../functions/includes/node/GetDofPosition.md)
- [GetVariable](../../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)