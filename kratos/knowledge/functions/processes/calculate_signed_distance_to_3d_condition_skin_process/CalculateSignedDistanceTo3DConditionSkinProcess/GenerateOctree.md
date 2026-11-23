---
type: C++ Method
title: GenerateOctree
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1225-L1297
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/SetBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateOctree()`

# Calls

- [Start](../../../../functions/utilities/timer/Timer/Start.md)
- [SetBoundingBox](../../../../functions/spatial_containers/octree_binary/OctreeBinary/SetBoundingBox.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)