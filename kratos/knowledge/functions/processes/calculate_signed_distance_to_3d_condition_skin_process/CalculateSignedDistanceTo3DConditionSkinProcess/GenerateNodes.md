---
type: C++ Method
title: GenerateNodes
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1307-L1328
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
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/pGetDataPointer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateNodes()`

# Calls

- [Start](../../../../functions/utilities/timer/Timer/Start.md)
- [GetAllLeavesVector](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetAllLeavesVector.md)
- [pGetDataPointer](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/pGetDataPointer.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)