---
type: C++ Method
title: PrintOutput
resource: input_output/vtk_output.cpp#L165-L190
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::PrintOutput(const std::string& rOutputFilename)`

# Calls

- [WriteModelPartToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GlobalNumberOfNodes](../../../../functions/sources/communicator/Communicator/GlobalNumberOfNodes.md)
- [GlobalNumberOfElements](../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)
- [WriteModelPartWithoutNodesToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile.md)