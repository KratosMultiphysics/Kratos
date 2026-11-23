---
type: C++ Method
title: GenerateStepLabel
resource: utilities/input_output_utilities.cpp#L92-L108
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string InputOutputUtilities::GenerateStepLabel( const std::size_t Step, const std::size_t Spaces )`

# Called by

- [PrintOutput](../../../../functions/input_output/ensight_output/EnSightOutput/PrintOutput.md)
- [WriteGeometryFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile.md)