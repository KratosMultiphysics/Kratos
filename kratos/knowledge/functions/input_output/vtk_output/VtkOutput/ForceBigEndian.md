---
type: C++ Method
title: ForceBigEndian
resource: input_output/vtk_output.cpp#L1033-L1043
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtk_output/WriteScalarDataToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/WriteVectorDataToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::ForceBigEndian(unsigned char* pBytes) const`

# Called by

- [WriteScalarDataToFile](../../../../functions/input_output/vtk_output/WriteScalarDataToFile.md)
- [WriteVectorDataToFile](../../../../functions/input_output/vtk_output/WriteVectorDataToFile.md)