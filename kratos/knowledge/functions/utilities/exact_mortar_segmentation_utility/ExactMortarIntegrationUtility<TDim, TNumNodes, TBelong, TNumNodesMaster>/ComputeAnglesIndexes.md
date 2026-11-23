---
type: C++ Method
title: ComputeAnglesIndexes
resource: utilities/exact_mortar_segmentation_utility.cpp#L1376-L1404
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/AnglePoints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TriangleIntersections
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline std::vector<IndexType> ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::ComputeAnglesIndexes( PointListType& rPointList, const array_1d<double, 3>& rNormal ) const`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [AnglePoints](../../../../functions/utilities/exact_mortar_segmentation_utility/AnglePoints.md)

# Called by

- [TriangleIntersections](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TriangleIntersections.md)