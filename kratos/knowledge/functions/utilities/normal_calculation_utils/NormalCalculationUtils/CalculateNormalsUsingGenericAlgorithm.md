---
type: C++ Method
title: CalculateNormalsUsingGenericAlgorithm
resource: utilities/normal_calculation_utils.cpp#L580-L636
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm( ModelPart& rModelPart, const bool ConsiderUnitNormal, const NormalVariableType& rNormalVariable )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)