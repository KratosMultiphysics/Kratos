---
type: C++ Method
title: WriteModelPartToFile
resource: input_output/vtk_output.cpp#L195-L219
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/WriteHeaderToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteMeshToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteModelPartToFile(const ModelPart& rModelPart, const bool IsSubModelPart, const std::string& rOutputFilename)`

# Calls

- [WriteHeaderToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteHeaderToFile.md)
- [WriteMeshToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteMeshToFile.md)
- [WriteNodalResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalResultsToFile.md)
- [WriteElementResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile.md)
- [WriteConditionResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile.md)

# Called by

- [PrintOutput](../../../../functions/input_output/vtk_output/VtkOutput/PrintOutput.md)
- [WriteModelPartWithoutNodesToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile.md)