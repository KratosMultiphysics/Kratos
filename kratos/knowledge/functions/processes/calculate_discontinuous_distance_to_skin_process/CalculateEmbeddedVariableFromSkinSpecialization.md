---
type: C++ Function
title: CalculateEmbeddedVariableFromSkinSpecialization
resource: processes/calculate_discontinuous_distance_to_skin_process.h#L523-L591
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetNonHistoricalVariableToZero
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEmbeddedVariableFromSkinSpecialization( const Variable<TVarType> &rVariable, const Variable<TVarType> &rEmbeddedVariable)`

# Calls

- [SolutionStepsDataHas](../../../functions/includes/node/SolutionStepsDataHas.md)
- [SetNonHistoricalVariableToZero](../../../functions/utilities/variable_utils/SetNonHistoricalVariableToZero.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)