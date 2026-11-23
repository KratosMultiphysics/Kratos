---
type: C++ Method
title: WriteScalarSolutionStepVariable
resource: input_output/vtk_output.cpp#L834-L847
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/WriteScalarDataToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalScalarValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteScalarSolutionStepVariable( const TContainerType& rContainer, const TVarType& rVariable, std::ofstream& rFileStream) const`

# Calls

- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [WriteScalarDataToFile](../../../../functions/input_output/vtk_output/WriteScalarDataToFile.md)

# Called by

- [WriteNodalScalarValues](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalScalarValues.md)