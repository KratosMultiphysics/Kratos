---
type: C++ Function
title: pAddDof
resource: includes/node.h#L1202-L1221
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/SortDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/SetNodalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/SetReaction
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable)`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [SortDofs](../../../functions/includes/node/SortDofs.md)
- [SetNodalData](../../../functions/includes/dof/Dof/SetNodalData.md)
- [SetReaction](../../../functions/includes/dof/SetReaction.md)