---
type: C++ Function
title: ApplyFixity
resource: utilities/variable_utils.h#L1088-L1120
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/CheckVariableExists
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FixDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FreeDof
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyFixity( const TVarType& rVar, const bool IsFixed, NodesContainerType& rNodes )`

# Calls

- [HasDofFor](../../../functions/includes/node/HasDofFor.md)
- [CheckVariableExists](../../../functions/utilities/variable_utils/CheckVariableExists.md)
- [FixDof](../../../functions/includes/dof/Dof/FixDof.md)
- [FreeDof](../../../functions/includes/dof/Dof/FreeDof.md)