---
type: C++ Method
title: ComputeUnitNormalsFromAreaNormals
resource: utilities/normal_calculation_utils.cpp#L365-L423
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::ComputeUnitNormalsFromAreaNormals( ModelPart& rModelPart, const Variable<array_1d<double,3>>& rNormalVariable)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)