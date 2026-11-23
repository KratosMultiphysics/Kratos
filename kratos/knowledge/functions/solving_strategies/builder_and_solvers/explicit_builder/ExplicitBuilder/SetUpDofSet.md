---
type: C++ Method
title: SetUpDofSet
resource: solving_strategies/builder_and_solvers/explicit_builder.h#L565-L658
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/HasReaction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void SetUpDofSet(const ModelPart& rModelPart)`

# Calls

- [HasReaction](../../../../../functions/includes/dof/Dof/HasReaction.md)
- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)