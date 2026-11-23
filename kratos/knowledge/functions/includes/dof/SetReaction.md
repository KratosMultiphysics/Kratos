---
type: C++ Function
title: SetReaction
resource: includes/dof.h#L298-L302
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/variables_list/SetDofReaction
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/node/pAddDof
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetReaction(TReactionType const& rReaction)`

# Calls

- [SetDofReaction](../../../functions/containers/variables_list/SetDofReaction.md)

# Called by

- [pAddDof](../../../functions/includes/node/pAddDof.md)