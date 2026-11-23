---
type: C++ Method
title: InitializeNormals
resource: utilities/normal_calculation_utils.cpp#L63-L105
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::InitializeNormals( ModelPart& rModelPart, const Variable<array_1d<double,3>>& rNormalVariable )`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)