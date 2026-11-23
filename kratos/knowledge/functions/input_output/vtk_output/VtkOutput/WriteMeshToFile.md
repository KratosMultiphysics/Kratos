---
type: C++ Method
title: WriteMeshToFile
resource: input_output/vtk_output.cpp#L316-L320
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodesToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteMeshToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const`

# Calls

- [WriteNodesToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodesToFile.md)
- [WriteConditionsAndElementsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile.md)

# Called by

- [WriteModelPartToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile.md)