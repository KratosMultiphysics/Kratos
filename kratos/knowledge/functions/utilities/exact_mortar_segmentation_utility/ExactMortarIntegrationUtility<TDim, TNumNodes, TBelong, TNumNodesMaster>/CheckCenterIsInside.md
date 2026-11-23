---
type: C++ Method
title: CheckCenterIsInside
resource: utilities/exact_mortar_segmentation_utility.cpp#L1587-L1609
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TriangleIntersections
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::CheckCenterIsInside( const array_1d<double, 2>& rAuxiliarCenterLocalCoordinates, const SizeType NumNodes )`

# Called by

- [TriangleIntersections](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TriangleIntersections.md)