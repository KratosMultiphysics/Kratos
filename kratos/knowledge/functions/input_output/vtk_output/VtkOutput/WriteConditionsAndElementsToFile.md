---
type: C++ Method
title: WriteConditionsAndElementsToFile
resource: input_output/vtk_output.cpp#L350-L374
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkCellListSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteConnectivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteCellType
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteMeshToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteConditionsAndElementsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const`

# Calls

- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [DetermineVtkCellListSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkCellListSize.md)
- [WriteConnectivity](../../../../functions/input_output/vtk_output/VtkOutput/WriteConnectivity.md)
- [WriteCellType](../../../../functions/input_output/vtk_output/VtkOutput/WriteCellType.md)

# Called by

- [WriteMeshToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteMeshToFile.md)