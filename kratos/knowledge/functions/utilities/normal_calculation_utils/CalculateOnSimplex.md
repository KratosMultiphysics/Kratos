---
type: C++ Function
title: CalculateOnSimplex
resource: utilities/normal_calculation_utils.h#L249-L297
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormal2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateOnSimplex( ModelPart& rModelPart, const std::size_t Dimension, const Variable<TValueType>& rVariable, const TValueType Zero, const NormalVariableType& rNormalVariable = NORMAL )`

# Calls

- [CalculateNormal2D](../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormal2D.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)