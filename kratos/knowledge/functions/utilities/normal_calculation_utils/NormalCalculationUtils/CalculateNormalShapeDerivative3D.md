---
type: C++ Method
title: CalculateNormalShapeDerivative3D
resource: utilities/normal_calculation_utils.cpp#L504-L554
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivativesOnSimplex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::CalculateNormalShapeDerivative3D( ConditionType& rCondition )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CalculateNormalShapeDerivativesOnSimplex](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivativesOnSimplex.md)