---
type: C++ Method
title: WriteGeometricalContainerResults
resource: input_output/vtk_output.cpp#L721-L757
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteGeometricalContainerResults( const std::string& rVariableName, const TContainerType& rContainer, std::ofstream& rFileStream) const`

# Calls

- [WriteScalarContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable.md)
- [WriteVectorContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable.md)

# Called by

- [WriteElementResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile.md)
- [WriteConditionResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile.md)