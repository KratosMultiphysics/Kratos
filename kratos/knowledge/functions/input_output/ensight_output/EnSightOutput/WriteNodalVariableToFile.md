---
type: C++ Method
title: WriteNodalVariableToFile
resource: input_output/ensight_output.cpp#L776-L1306
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

`void EnSightOutput::WriteNodalVariableToFile( const std::string& rVariableName, const bool IsHistorical, const std::string& rFileName )`

# Calls

- [GetExtensionFile](../../../../functions/input_output/ensight_output/EnSightOutput/GetExtensionFile.md)
- [WriteString](../../../../functions/input_output/ensight_output/EnSightOutput/WriteString.md)
- [GetTypeLabel](../../../../functions/input_output/ensight_output/EnSightOutput/GetTypeLabel.md)
- [WriteScalarData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteScalarData.md)
- [WriteVectorData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteVectorData.md)
- [WriteSymmetricTensorData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteSymmetricTensorData.md)
- [WriteAsymmetricTensorData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteAsymmetricTensorData.md)

# Called by

- [PrintOutput](../../../../functions/input_output/ensight_output/EnSightOutput/PrintOutput.md)