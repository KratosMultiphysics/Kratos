---
type: C++ Function
title: ProjectIterativeLine2D
resource: utilities/geometrical_projection_utilities.h#L319-L395
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<2, 2, false>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<2, 2, true>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool ProjectIterativeLine2D( TGeometryType& rGeomOrigin, const array_1d<double,3>& rPointDestiny, array_1d<double,3>& rResultingPoint, const array_1d<double, 3>& rNormal, const double Tolerance = 1.0e-8, double DeltaXi = 0.5 )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<2, 2, false>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<2, 2, true>/GetExactIntegration.md)