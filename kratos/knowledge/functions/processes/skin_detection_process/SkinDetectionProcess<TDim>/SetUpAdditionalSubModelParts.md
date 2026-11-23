---
type: C++ Method
title: SetUpAdditionalSubModelParts
resource: processes/skin_detection_process.cpp#L321-L397
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SkinDetectionProcess<TDim>::SetUpAdditionalSubModelParts(const ModelPart& rAuxiliaryModelPart)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [Execute](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute.md)