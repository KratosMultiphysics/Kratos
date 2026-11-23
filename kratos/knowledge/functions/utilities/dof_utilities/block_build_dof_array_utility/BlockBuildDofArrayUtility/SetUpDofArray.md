---
type: C++ Method
title: SetUpDofArray
resource: utilities/dof_utilities/block_build_dof_array_utility.h#L69-L178
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/HasReaction
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void SetUpDofArray( const ModelPart& rModelPart, DofsArrayType& rDofArray, const unsigned int EchoLevel = 0, const bool CheckReactionDofs = false)`

# Calls

- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)
- [FullName](../../../../../functions/includes/model_part/FullName.md)
- [HasReaction](../../../../../functions/includes/dof/Dof/HasReaction.md)