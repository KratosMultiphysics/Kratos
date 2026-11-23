---
type: C++ Method
title: WritePropertiesIdsToFile
resource: input_output/vtk_output.cpp#L1100-L1112
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
  - target: functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WritePropertiesIdsToFile( const TContainerType& rContainer, std::ofstream& rFileStream) const`

# Calls

- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [WriteScalarDataToFile](../../../../functions/input_output/vtk_output/WriteScalarDataToFile.md)

# Called by

- [WriteElementResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile.md)
- [WriteConditionResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile.md)