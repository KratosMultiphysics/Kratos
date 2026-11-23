---
type: C++ Method
title: WriteGeometricalGaussVariableToFile
resource: input_output/ensight_output.cpp#L1785-L2137
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/ensight_output/EnSightOutput/GetExtensionFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/GetTypeLabel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteScalarData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteVectorData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteSymmetricTensorData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteAsymmetricTensorData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EnSightOutput::WriteGeometricalGaussVariableToFile( const std::string& rVariableName, const std::string& rFileName, const bool IsElement )`

# Calls

- [GetExtensionFile](../../../../functions/input_output/ensight_output/EnSightOutput/GetExtensionFile.md)
- [WriteString](../../../../functions/input_output/ensight_output/EnSightOutput/WriteString.md)
- [GetTypeLabel](../../../../functions/input_output/ensight_output/EnSightOutput/GetTypeLabel.md)
- [WriteScalarData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteScalarData.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [WriteVectorData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteVectorData.md)
- [WriteSymmetricTensorData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteSymmetricTensorData.md)
- [WriteAsymmetricTensorData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteAsymmetricTensorData.md)

# Called by

- [PrintOutput](../../../../functions/input_output/ensight_output/EnSightOutput/PrintOutput.md)