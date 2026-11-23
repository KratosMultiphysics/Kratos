---
type: C++ Method
title: FindOptimumPositionsAndWeights
resource: processes/tetrahedra_mesh_worst_element_smoothing_process.cpp#L46-L65
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/tetrahedra_mesh_worst_element_smoothing_process/TetrahedraMeshWorstElementSmoothingProcess/CalculateElementOptimumPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/Quality
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TetrahedraMeshWorstElementSmoothingProcess::FindOptimumPositionsAndWeights(NodeType& rNode, PointsVectorType& rOptimumPoints, Vector& rWeights)`

# Calls

- [CalculateElementOptimumPosition](../../../../functions/processes/tetrahedra_mesh_worst_element_smoothing_process/TetrahedraMeshWorstElementSmoothingProcess/CalculateElementOptimumPosition.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Quality](../../../../functions/geometries/geometry/Quality.md)