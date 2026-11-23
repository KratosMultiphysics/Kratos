---
type: C++ Method
title: ApplyProxy
resource: utilities/pointer_map_communicator.h#L115-L133
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/pointer_map_communicator/GlobalPointerMapCommunicator/GetMyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/openmp_utils/IsInParallel
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ApplyProxy( const TApplyFunctor& rApplyFunctor, TGPMapCommunicator& rPointerCommunicator) : mCurrentRank(rPointerCommunicator.GetMyPID()), mrApplyFunctor(rApplyFunctor), mrPointerCommunicator(rPointerCommunicator)`

# Calls

- [GetMyPID](../../../../functions/utilities/pointer_map_communicator/GlobalPointerMapCommunicator/GetMyPID.md)
- [IsInParallel](../../../../functions/utilities/openmp_utils/IsInParallel.md)