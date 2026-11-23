---
type: C++ Method
title: SetUpDofArray
resource: utilities/dof_utilities/elimination_build_dof_array_utility.h#L77-L172
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/openmp_utils/ThisThread
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

`static void SetUpDofArray( const ModelPart& rModelPart, DofsArrayType& rDofArray, #ifdef USE_LOCKS_IN_ASSEMBLY std::vector<omp_lock_t>& rLockArray, #endif const unsigned int EchoLevel = 0, const bool CheckReactionDofs = false)`

# Calls

- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)
- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)
- [FullName](../../../../../functions/includes/model_part/FullName.md)
- [HasReaction](../../../../../functions/includes/dof/Dof/HasReaction.md)