---
type: C++ Method
title: GetExactIntegration
resource: utilities/exact_mortar_segmentation_utility.cpp#L29-L134
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_projection_utilities/FastProjectDirection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/ProjectIterativeLine2D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ExactMortarIntegrationUtility<2, 2, false>::GetExactIntegration( const GeometryType& rOriginalSlaveGeometry, const array_1d<double, 3>& rSlaveNormal, const GeometryType& rOriginalMasterGeometry, const array_1d<double, 3>& rMasterNormal, ConditionArrayListType& rConditionsPointsSlave )`

# Calls

- [FastProjectDirection](../../../../functions/utilities/geometrical_projection_utilities/FastProjectDirection.md)
- [ProjectIterativeLine2D](../../../../functions/utilities/geometrical_projection_utilities/ProjectIterativeLine2D.md)