---
type: C++ Function
title: CheckPoints
resource: utilities/exact_mortar_segmentation_utility.h#L535-L543
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeClippingIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/PushBackPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool CheckPoints( const Point& rPointOrig, const Point& rPointDest )`

# Called by

- [ComputeClippingIntersections](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeClippingIntersections.md)
- [PushBackPoints](../../../functions/utilities/exact_mortar_segmentation_utility/PushBackPoints.md)