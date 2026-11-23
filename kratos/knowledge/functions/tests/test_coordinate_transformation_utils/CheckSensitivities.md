---
type: Python Function
title: CheckSensitivities
resource: tests/test_coordinate_transformation_utils.py#L66-L77
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/coordinate_transformation_utilities/CoordinateTransformationUtils/CalculateRotationOperatorPure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_coordinate_transformation_utils/FiniteDifferenceNormalShapeSensitivityTest
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CheckSensitivities(node, derivative_node_index, derivative_direction_index, ref_value): # calculate analytical shape sensitivities`

# Calls

- [CalculateRotationOperatorPure](../../../functions/utilities/coordinate_transformation_utilities/CoordinateTransformationUtils/CalculateRotationOperatorPure.md)
- [assertMatrixAlmostEqual](../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)

# Called by

- [FiniteDifferenceNormalShapeSensitivityTest](../../../functions/tests/test_coordinate_transformation_utils/FiniteDifferenceNormalShapeSensitivityTest.md)