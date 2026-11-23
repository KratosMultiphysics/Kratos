---
type: C++ Method
title: WeightedAccumulateVariableOnNodes
resource: utilities/variable_utils.cpp#L344-L389
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

`void VariableUtils::WeightedAccumulateVariableOnNodes( ModelPart& rModelPart, const Variable<TDataType>& rVariable, const Variable<TWeightDataType>& rWeightVariable, const bool IsInverseWeightProvided)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)