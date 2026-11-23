---
type: C++ Method
title: SetUpDofArray
resource: utilities/dof_utilities/dof_array_utilities.cpp#L26-L100
visibility: private
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
---

# Signature

`void DofArrayUtilities::SetUpDofArray( const ModelPart& rModelPart, DofsArrayType& rDofArray, const unsigned int EchoLevel)`

# Calls

- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)
- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)
- [FullName](../../../../../functions/includes/model_part/FullName.md)