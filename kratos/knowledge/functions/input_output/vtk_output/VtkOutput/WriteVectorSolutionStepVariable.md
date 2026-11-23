---
type: C++ Method
title: WriteVectorSolutionStepVariable
resource: input_output/vtk_output.cpp#L853-L873
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/WriteVectorDataToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalVectorValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteVectorSolutionStepVariable( const TContainerType& rContainer, const TVarType& rVariable, std::ofstream& rFileStream) const`

# Calls

- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [WriteVectorDataToFile](../../../../functions/input_output/vtk_output/WriteVectorDataToFile.md)

# Called by

- [WriteNodalVectorValues](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalVectorValues.md)