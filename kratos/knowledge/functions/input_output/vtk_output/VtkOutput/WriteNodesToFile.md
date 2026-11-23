---
type: C++ Method
title: WriteNodesToFile
resource: input_output/vtk_output.cpp#L325-L345
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/WriteVectorDataToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteMeshToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteNodesToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const`

# Calls

- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [WriteVectorDataToFile](../../../../functions/input_output/vtk_output/WriteVectorDataToFile.md)

# Called by

- [WriteMeshToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteMeshToFile.md)