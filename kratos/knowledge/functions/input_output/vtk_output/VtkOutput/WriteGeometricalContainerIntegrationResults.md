---
type: C++ Method
title: WriteGeometricalContainerIntegrationResults
resource: input_output/vtk_output.cpp#L763-L790
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable
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

`void VtkOutput::WriteGeometricalContainerIntegrationResults( const std::string& rVariableName, const TContainerType& rContainer, std::ofstream& rFileStream) const`

# Calls

- [WriteIntegrationScalarContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable.md)
- [WriteIntegrationVectorContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable.md)

# Called by

- [WriteElementResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile.md)
- [WriteConditionResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile.md)