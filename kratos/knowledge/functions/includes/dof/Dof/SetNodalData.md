---
type: C++ Method
title: SetNodalData
resource: includes/dof.h#L351-L361
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list/pGetDofReaction
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/node/Node
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/pAddDof
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetNodalData(NodalData* pNewNodalData)`

# Calls

- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [pGetDofReaction](../../../../functions/containers/variables_list/pGetDofReaction.md)

# Called by

- [Node](../../../../functions/includes/node/Node.md)
- [pAddDof](../../../../functions/includes/node/pAddDof.md)