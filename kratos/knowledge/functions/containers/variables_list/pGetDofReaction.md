---
type: C++ Function
title: pGetDofReaction
resource: containers/variables_list.h#L327-L329
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/dof/Dof/SetNodalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/HasReaction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_variables_list/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`const VariableData* pGetDofReaction(int DofIndex) const`

# Called by

- [SetNodalData](../../../functions/includes/dof/Dof/SetNodalData.md)
- [HasReaction](../../../functions/includes/dof/Dof/HasReaction.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/containers/test_variables_list/KRATOS_TEST_CASE_IN_SUITE.md)