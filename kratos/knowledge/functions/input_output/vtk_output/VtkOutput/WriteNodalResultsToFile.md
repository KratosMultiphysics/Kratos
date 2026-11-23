---
type: C++ Method
title: WriteNodalResultsToFile
resource: input_output/vtk_output.cpp#L490-L554
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/IsCompatibleVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalContainerResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIdsToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteNodalResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)`

# Calls

- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [IsCompatibleVariable](../../../../functions/input_output/vtk_output/VtkOutput/IsCompatibleVariable.md)
- [WriteNodalContainerResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalContainerResults.md)
- [WriteFlagContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable.md)
- [WriteIdsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteIdsToFile.md)

# Called by

- [WriteModelPartToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile.md)