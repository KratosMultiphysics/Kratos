---
type: C++ Method
title: CalculateRotationOperatorPure
resource: utilities/coordinate_transformation_utilities.h#L119-L142
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_coordinate_transformation_utils/FiniteDifferenceNormalShapeSensitivityTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_coordinate_transformation_utils/CheckSensitivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateRotationOperatorPure( TLocalMatrixType& rRotationMatrix, const GeometryType::PointType& rThisPoint) const`

# Called by

- [FiniteDifferenceNormalShapeSensitivityTest](../../../../functions/tests/test_coordinate_transformation_utils/FiniteDifferenceNormalShapeSensitivityTest.md)
- [CheckSensitivities](../../../../functions/tests/test_coordinate_transformation_utils/CheckSensitivities.md)