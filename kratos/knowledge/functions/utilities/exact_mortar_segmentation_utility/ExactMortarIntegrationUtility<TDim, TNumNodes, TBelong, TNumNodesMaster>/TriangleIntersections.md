---
type: C++ Method
title: TriangleIntersections
resource: utilities/exact_mortar_segmentation_utility.cpp#L1459-L1581
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeClippingIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/RotatePoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeAnglesIndexes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/CheckCenterIsInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::TriangleIntersections( ConditionArrayListType& rConditionsPointsSlave, PointListType& rPointList, const TGeometryType& rOriginalSlaveGeometry, const GeometryPointType& rSlaveGeometry, const GeometryPointType& rMasterGeometry, const array_1d<double, 3>& rSlaveTangentXi, const array_1d<double, 3>& rSlaveTangentEta, const Point& rRefCenter, const bool IsAllInside )`

# Calls

- [ComputeClippingIntersections](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeClippingIntersections.md)
- [RotatePoint](../../../../functions/utilities/mortar_utilities/RotatePoint.md)
- [ComputeAnglesIndexes](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeAnglesIndexes.md)
- [CheckCenterIsInside](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/CheckCenterIsInside.md)