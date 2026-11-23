---
type: C++ Method
title: PerformSmoothing
resource: processes/mesh_local_smoothing_process.cpp#L122-L136
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/InterpolateNodeOptimumPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/MoveNodeIfImprovesMinimumQuality
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshLocalSmoothingProcess::PerformSmoothing()`

# Calls

- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [InterpolateNodeOptimumPosition](../../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/InterpolateNodeOptimumPosition.md)
- [MoveNodeIfImprovesMinimumQuality](../../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/MoveNodeIfImprovesMinimumQuality.md)

# Called by

- [Execute](../../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/Execute.md)