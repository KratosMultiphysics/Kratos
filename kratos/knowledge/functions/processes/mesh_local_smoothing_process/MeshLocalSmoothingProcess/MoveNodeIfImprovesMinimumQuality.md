---
type: C++ Method
title: MoveNodeIfImprovesMinimumQuality
resource: processes/mesh_local_smoothing_process.cpp#L152-L202
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/tetrahedra_ball/TetrahedraBall/CalculateMinQuality
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PerformSmoothing
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshLocalSmoothingProcess::MoveNodeIfImprovesMinimumQuality(NodeType& rNode, Point const& OptimumPosition)`

# Calls

- [CalculateMinQuality](../../../../functions/modeler/tetrahedra_ball/TetrahedraBall/CalculateMinQuality.md)

# Called by

- [PerformSmoothing](../../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PerformSmoothing.md)