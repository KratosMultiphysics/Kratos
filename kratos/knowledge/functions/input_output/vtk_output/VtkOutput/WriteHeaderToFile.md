---
type: C++ Method
title: WriteHeaderToFile
resource: input_output/vtk_output.cpp#L298-L311
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteHeaderToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const`

# Called by

- [WriteModelPartToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile.md)