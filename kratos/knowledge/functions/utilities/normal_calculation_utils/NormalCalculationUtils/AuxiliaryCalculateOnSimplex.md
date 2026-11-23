---
type: C++ Method
title: AuxiliaryCalculateOnSimplex
resource: utilities/normal_calculation_utils.cpp#L726-L777
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
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormal2D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::AuxiliaryCalculateOnSimplex( ConditionsArrayType& rConditions, const std::size_t Dimension, const NormalVariableType& rNormalVariable )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [CalculateNormal2D](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormal2D.md)