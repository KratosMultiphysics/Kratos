---
type: C++ Function
title: Quality
resource: geometries/geometry.h#L1766-L1798
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/tetrahedra_ball/TetrahedraBall/CalculateMinQuality
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/SelectLowQualityElementNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedra_mesh_quality_weighted_smoothing_process/TetrahedraMeshQualityWeightedSmoothingProcess/FindOptimumPositionsAndWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedra_mesh_worst_element_smoothing_process/TetrahedraMeshWorstElementSmoothingProcess/FindOptimumPositionsAndWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_tetrahedra_3d_4/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double Quality(const QualityCriteria qualityCriteria) const`

# Called by

- [CalculateMinQuality](../../../functions/modeler/tetrahedra_ball/TetrahedraBall/CalculateMinQuality.md)
- [SelectLowQualityElementNodes](../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/SelectLowQualityElementNodes.md)
- [FindOptimumPositionsAndWeights](../../../functions/processes/tetrahedra_mesh_quality_weighted_smoothing_process/TetrahedraMeshQualityWeightedSmoothingProcess/FindOptimumPositionsAndWeights.md)
- [FindOptimumPositionsAndWeights](../../../functions/processes/tetrahedra_mesh_worst_element_smoothing_process/TetrahedraMeshWorstElementSmoothingProcess/FindOptimumPositionsAndWeights.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_tetrahedra_3d_4/KRATOS_TEST_CASE_IN_SUITE.md)