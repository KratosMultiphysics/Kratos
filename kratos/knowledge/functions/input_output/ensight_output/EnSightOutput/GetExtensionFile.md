---
type: C++ Method
title: GetExtensionFile
resource: input_output/ensight_output.cpp#L2745-L2764
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/WriteCaseFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteNodalVariableToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteNodalFlagToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometricalVariableToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometricalFlagToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometricalGaussVariableToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string EnSightOutput::GetExtensionFile(const VariableType Type) const`

# Called by

- [WriteCaseFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteCaseFile.md)
- [WriteNodalVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteNodalVariableToFile.md)
- [WriteNodalFlagToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteNodalFlagToFile.md)
- [WriteGeometricalVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalVariableToFile.md)
- [WriteGeometricalFlagToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalFlagToFile.md)
- [WriteGeometricalGaussVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalGaussVariableToFile.md)