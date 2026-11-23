---
type: C++ Function
title: CheckInside
resource: utilities/exact_mortar_segmentation_utility.h#L638-L649
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false, 4>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false, 3>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true, 4>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true, 3>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void CheckInside( array_1d<bool, TSizeCheck>& rAllInside, GeometryPointType& rGeometry1, GeometryPointType& rGeometry2, const double Tolerance )`

# Called by

- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false, 4>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false, 3>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true, 4>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true, 3>/GetExactIntegration.md)