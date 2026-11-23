---
type: C++ Method
title: GetExactIntegration
resource: utilities/exact_mortar_segmentation_utility.cpp#L225-L301
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastProject
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/RotatePoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/CheckInside
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/CheckAllInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ExactMortarIntegrationUtility<3, 4, false>::GetExactIntegration( const GeometryType& rOriginalSlaveGeometry, const array_1d<double, 3>& rSlaveNormal, const GeometryType& rOriginalMasterGeometry, const array_1d<double, 3>& rMasterNormal, ConditionArrayListType& rConditionsPointsSlave )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [FastProject](../../../../functions/utilities/geometrical_projection_utilities/FastProject.md)
- [RotatePoint](../../../../functions/utilities/mortar_utilities/RotatePoint.md)
- [CheckInside](../../../../functions/utilities/exact_mortar_segmentation_utility/CheckInside.md)
- [CheckAllInside](../../../../functions/utilities/exact_mortar_segmentation_utility/CheckAllInside.md)