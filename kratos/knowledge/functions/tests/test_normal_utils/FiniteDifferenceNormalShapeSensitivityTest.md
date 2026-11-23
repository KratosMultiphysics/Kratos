---
type: Python Function
title: FiniteDifferenceNormalShapeSensitivityTest
resource: tests/test_normal_utils.py#L49-L95
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivativesOnSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/check_normal_sensitivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def FiniteDifferenceNormalShapeSensitivityTest(UnitTestObject, model_part, check_condition_ids, delta, tolerance):`

# Calls

- [CalculateNormalShapeDerivativesOnSimplex](../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivativesOnSimplex.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [check_normal_sensitivities](../../../functions/tests/test_normal_utils/check_normal_sensitivities.md)