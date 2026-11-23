---
type: C++ Method
title: EvaluateNumberOfSubsteps
resource: processes/levelset_convection_process.h#L622-L668
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`unsigned int EvaluateNumberOfSubsteps()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [Execute](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute.md)