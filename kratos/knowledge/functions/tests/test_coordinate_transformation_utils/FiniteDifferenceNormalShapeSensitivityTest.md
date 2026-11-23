---
type: Python Function
title: FiniteDifferenceNormalShapeSensitivityTest
resource: tests/test_coordinate_transformation_utils.py#L36-L98
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivativesOnSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/CoordinateTransformationUtils/CalculateRotationOperatorPure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_coordinate_transformation_utils/CheckSensitivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def FiniteDifferenceNormalShapeSensitivityTest(UnitTestObject, model_part, check_node_ids, delta, tolerance):`

# Calls

- [CalculateNormalShapeDerivativesOnSimplex](../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivativesOnSimplex.md)
- [CalculateRotationOperatorPure](../../../functions/utilities/coordinate_transformation_utilities/CoordinateTransformationUtils/CalculateRotationOperatorPure.md)
- [CheckSensitivities](../../../functions/tests/test_coordinate_transformation_utils/CheckSensitivities.md)