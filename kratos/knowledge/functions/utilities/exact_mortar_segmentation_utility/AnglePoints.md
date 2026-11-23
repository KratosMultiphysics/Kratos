---
type: C++ Function
title: AnglePoints
resource: utilities/exact_mortar_segmentation_utility.h#L512-L527
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeAnglesIndexes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double AnglePoints( const Point& rPointOrig1, const Point& rPointOrig2, const array_1d<double, 3>& rAxis1, const array_1d<double, 3>& rAxis2 )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [ComputeAnglesIndexes](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeAnglesIndexes.md)