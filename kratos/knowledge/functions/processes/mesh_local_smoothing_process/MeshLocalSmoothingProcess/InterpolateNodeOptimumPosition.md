---
type: C++ Method
title: InterpolateNodeOptimumPosition
resource: processes/mesh_local_smoothing_process.cpp#L204-L221
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PerformSmoothing
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshLocalSmoothingProcess::InterpolateNodeOptimumPosition(PointsVectorType const& rOptimumPoints, Vector const& rWeights, Point& OptimumPosition)`

# Called by

- [PerformSmoothing](../../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PerformSmoothing.md)