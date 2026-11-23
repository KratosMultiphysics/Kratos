---
type: C++ Method
title: CalculateElementOptimumPosition
resource: processes/tetrahedra_mesh_worst_element_smoothing_process.cpp#L67-L83
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/tetrahedra_mesh_quality_weighted_smoothing_process/TetrahedraMeshQualityWeightedSmoothingProcess/FindOptimumPositionsAndWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedra_mesh_worst_element_smoothing_process/TetrahedraMeshWorstElementSmoothingProcess/FindOptimumPositionsAndWeights
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TetrahedraMeshWorstElementSmoothingProcess::CalculateElementOptimumPosition(NodeType& rNode, Geometry<Node > const& rTetrahedra, Point& rOptimumPoint)`

# Called by

- [FindOptimumPositionsAndWeights](../../../../functions/processes/tetrahedra_mesh_quality_weighted_smoothing_process/TetrahedraMeshQualityWeightedSmoothingProcess/FindOptimumPositionsAndWeights.md)
- [FindOptimumPositionsAndWeights](../../../../functions/processes/tetrahedra_mesh_worst_element_smoothing_process/TetrahedraMeshWorstElementSmoothingProcess/FindOptimumPositionsAndWeights.md)