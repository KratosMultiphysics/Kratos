---
type: C++ Method
title: WriteConditionResultsToFile
resource: input_output/vtk_output.cpp#L620-L676
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/IsCompatibleVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WritePropertiesIdsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIdsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerIntegrationResults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteConditionResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [IsCompatibleVariable](../../../../functions/input_output/vtk_output/VtkOutput/IsCompatibleVariable.md)
- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [WriteGeometricalContainerResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerResults.md)
- [WriteFlagContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable.md)
- [WritePropertiesIdsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WritePropertiesIdsToFile.md)
- [WriteIdsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteIdsToFile.md)
- [WriteGeometricalContainerIntegrationResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerIntegrationResults.md)

# Called by

- [WriteModelPartToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile.md)