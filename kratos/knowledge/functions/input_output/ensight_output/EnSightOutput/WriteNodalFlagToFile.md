---
type: C++ Method
title: WriteNodalFlagToFile
resource: input_output/ensight_output.cpp#L1311-L1383
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
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EnSightOutput::WriteNodalFlagToFile( const std::string& rFlagName, const std::string& rFileName )`

# Calls

- [GetExtensionFile](../../../../functions/input_output/ensight_output/EnSightOutput/GetExtensionFile.md)
- [WriteString](../../../../functions/input_output/ensight_output/EnSightOutput/WriteString.md)
- [GetTypeLabel](../../../../functions/input_output/ensight_output/EnSightOutput/GetTypeLabel.md)
- [WriteScalarData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteScalarData.md)
- [IsDefined](../../../../functions/containers/flags/IsDefined.md)

# Called by

- [PrintOutput](../../../../functions/input_output/ensight_output/EnSightOutput/PrintOutput.md)