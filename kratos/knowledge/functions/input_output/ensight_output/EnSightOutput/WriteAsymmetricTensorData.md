---
type: C++ Method
title: WriteAsymmetricTensorData
resource: input_output/ensight_output.cpp#L2380-L2422
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/WriteNodalVariableToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometricalVariableToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometricalGaussVariableToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EnSightOutput::WriteAsymmetricTensorData( std::ofstream& rFileStream, const TData& rData, const bool EndOfLine, const bool AddInitialTabulation, const bool AddEndTabulation ) const`

# Called by

- [WriteNodalVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteNodalVariableToFile.md)
- [WriteGeometricalVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalVariableToFile.md)
- [WriteGeometricalGaussVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalGaussVariableToFile.md)