---
type: C++ Method
title: WriteString
resource: input_output/ensight_output.cpp#L2430-L2449
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile
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

`void EnSightOutput::WriteString( std::ofstream& rFileStream, const std::string& rString, const bool EndOfLine, const bool AddInitialTabulation, const bool AddEndTabulation ) const`

# Called by

- [WriteGeometryFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile.md)
- [WriteNodalVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteNodalVariableToFile.md)
- [WriteNodalFlagToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteNodalFlagToFile.md)
- [WriteGeometricalVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalVariableToFile.md)
- [WriteGeometricalFlagToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalFlagToFile.md)
- [WriteGeometricalGaussVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalGaussVariableToFile.md)