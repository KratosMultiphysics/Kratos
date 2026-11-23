---
type: C++ Method
title: CalculateNormal3D
resource: utilities/normal_calculation_utils.cpp#L477-L499
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::CalculateNormal3D( Condition& rCondition, array_1d<double,3>& rAn, array_1d<double,3>& rv1, array_1d<double,3>& rv2, const Variable<array_1d<double,3>>& rNormalVariable )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)