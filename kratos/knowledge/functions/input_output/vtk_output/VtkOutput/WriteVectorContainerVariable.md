---
type: C++ Method
title: WriteVectorContainerVariable
resource: input_output/vtk_output.cpp#L960-L980
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
  - target: functions/input_output/vtk_output/WriteVectorDataToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalVectorValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteVectorContainerVariable( const TContainerType& rContainer, const TVarType& rVariable, std::ofstream& rFileStream) const`

# Calls

- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [WriteVectorDataToFile](../../../../functions/input_output/vtk_output/WriteVectorDataToFile.md)

# Called by

- [WriteGeometricalContainerResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerResults.md)
- [WriteNodalVectorValues](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalVectorValues.md)