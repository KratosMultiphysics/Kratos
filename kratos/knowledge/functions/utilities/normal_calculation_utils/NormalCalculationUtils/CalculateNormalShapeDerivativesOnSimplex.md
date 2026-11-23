---
type: C++ Method
title: CalculateNormalShapeDerivativesOnSimplex
resource: utilities/normal_calculation_utils.cpp#L230-L261
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivative2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivative3D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_coordinate_transformation_utils/FiniteDifferenceNormalShapeSensitivityTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/FiniteDifferenceNormalShapeSensitivityTest
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::CalculateNormalShapeDerivativesOnSimplex( ConditionsArrayType& rConditions, const std::size_t Dimension )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CalculateNormalShapeDerivative2D](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivative2D.md)
- [CalculateNormalShapeDerivative3D](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivative3D.md)

# Called by

- [FiniteDifferenceNormalShapeSensitivityTest](../../../../functions/tests/test_coordinate_transformation_utils/FiniteDifferenceNormalShapeSensitivityTest.md)
- [FiniteDifferenceNormalShapeSensitivityTest](../../../../functions/tests/test_normal_utils/FiniteDifferenceNormalShapeSensitivityTest.md)