---
type: C++ Function
title: SortDofs
resource: includes/node.h#L1527-L1532
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/node/pAddDof
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SortDofs()`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)

# Called by

- [pAddDof](../../../functions/includes/node/pAddDof.md)