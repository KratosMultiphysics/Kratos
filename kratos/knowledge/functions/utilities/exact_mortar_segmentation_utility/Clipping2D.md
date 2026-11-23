---
type: C++ Function
title: Clipping2D
resource: utilities/exact_mortar_segmentation_utility.h#L445-L486
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeClippingIntersections
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool Clipping2D( Point& rPointIntersection, const Point& rPointOrig1, const Point& rPointOrig2, const Point& rPointDest1, const Point& rPointDest2 )`

# Called by

- [ComputeClippingIntersections](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeClippingIntersections.md)