---
type: C++ Method
title: ComputeClippingIntersections
resource: utilities/exact_mortar_segmentation_utility.cpp#L1410-L1452
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/exact_mortar_segmentation_utility/Clipping2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/CheckPoints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TriangleIntersections
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::ComputeClippingIntersections( PointListType& rPointList, const GeometryPointType& rSlaveGeometry, const GeometryPointType& rMasterGeometry, const Point& rRefCenter )`

# Calls

- [Clipping2D](../../../../functions/utilities/exact_mortar_segmentation_utility/Clipping2D.md)
- [CheckPoints](../../../../functions/utilities/exact_mortar_segmentation_utility/CheckPoints.md)

# Called by

- [TriangleIntersections](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TriangleIntersections.md)