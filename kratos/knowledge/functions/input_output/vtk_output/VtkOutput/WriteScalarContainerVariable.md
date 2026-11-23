---
type: C++ Method
title: WriteScalarContainerVariable
resource: input_output/vtk_output.cpp#L900-L914
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/WriteScalarDataToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalScalarValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteScalarContainerVariable( const TContainerType& rContainer, const TVarType& rVariable, std::ofstream& rFileStream) const`

# Calls

- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [WriteScalarDataToFile](../../../../functions/input_output/vtk_output/WriteScalarDataToFile.md)

# Called by

- [WriteGeometricalContainerResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerResults.md)
- [WriteNodalScalarValues](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalScalarValues.md)