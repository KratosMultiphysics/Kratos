---
type: C++ Method
title: WriteGeometryFile
resource: input_output/ensight_output.cpp#L561-L771
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/input_output_utilities/InputOutputUtilities/GenerateStepLabel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/ComputeBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteScalarData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/GetGeometryConnectivity
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/EnSightOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EnSightOutput::WriteGeometryFile(const std::string& rFileName)`

# Calls

- [GenerateStepLabel](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/GenerateStepLabel.md)
- [WriteString](../../../../functions/input_output/ensight_output/EnSightOutput/WriteString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [ComputeBoundingBox](../../../../functions/input_output/ensight_output/EnSightOutput/ComputeBoundingBox.md)
- [WriteScalarData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteScalarData.md)
- [GetGeometryConnectivity](../../../../functions/input_output/ensight_output/EnSightOutput/GetGeometryConnectivity.md)

# Called by

- [EnSightOutput](../../../../functions/input_output/ensight_output/EnSightOutput/EnSightOutput.md)
- [PrintOutput](../../../../functions/input_output/ensight_output/EnSightOutput/PrintOutput.md)