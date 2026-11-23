---
type: C++ Method
title: CalculateDistance
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1915-L1947
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetAllLeavesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistance()`

# Calls

- [Start](../../../../functions/utilities/timer/Timer/Start.md)
- [GetAllLeavesVector](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetAllLeavesVector.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)