---
type: C++ Method
title: CalculateMinQuality
resource: modeler/tetrahedra_ball.cpp#L39-L49
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/Quality
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/MoveNodeIfImprovesMinimumQuality
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double TetrahedraBall::CalculateMinQuality(const Geometry<NodeType>::QualityCriteria QualityCriteria) const`

# Calls

- [Quality](../../../../functions/geometries/geometry/Quality.md)

# Called by

- [MoveNodeIfImprovesMinimumQuality](../../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/MoveNodeIfImprovesMinimumQuality.md)