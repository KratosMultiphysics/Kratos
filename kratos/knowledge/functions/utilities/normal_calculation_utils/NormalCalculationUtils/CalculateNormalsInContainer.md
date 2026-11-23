---
type: C++ Method
title: CalculateNormalsInContainer
resource: utilities/normal_calculation_utils.cpp#L27-L57
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

`void NormalCalculationUtils::CalculateNormalsInContainer( ModelPart& rModelPart, const Variable<array_1d<double,3>>& rNormalVariable )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)