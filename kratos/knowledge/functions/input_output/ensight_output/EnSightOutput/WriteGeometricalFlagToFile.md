---
type: C++ Method
title: WriteGeometricalFlagToFile
resource: input_output/ensight_output.cpp#L1691-L1780
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

`void EnSightOutput::WriteGeometricalFlagToFile( const std::string& rFlagName, const std::string& rFileName, const bool IsElement )`

# Calls

- [GetExtensionFile](../../../../functions/input_output/ensight_output/EnSightOutput/GetExtensionFile.md)
- [WriteString](../../../../functions/input_output/ensight_output/EnSightOutput/WriteString.md)
- [GetTypeLabel](../../../../functions/input_output/ensight_output/EnSightOutput/GetTypeLabel.md)
- [WriteScalarData](../../../../functions/input_output/ensight_output/EnSightOutput/WriteScalarData.md)
- [IsDefined](../../../../functions/containers/flags/IsDefined.md)

# Called by

- [PrintOutput](../../../../functions/input_output/ensight_output/EnSightOutput/PrintOutput.md)