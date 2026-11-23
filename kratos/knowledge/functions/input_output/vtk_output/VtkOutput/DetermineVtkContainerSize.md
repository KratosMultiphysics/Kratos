---
type: C++ Method
title: DetermineVtkContainerSize
resource: input_output/vtk_output.cpp#L380-L389
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteScalarSolutionStepVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteVectorSolutionStepVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WritePropertiesIdsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIdsToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t VtkOutput::DetermineVtkContainerSize(const TContainerType& rContainer) const`

# Calls

- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)

# Called by

- [WriteConditionsAndElementsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile.md)
- [WriteElementResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile.md)
- [WriteConditionResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile.md)
- [WriteScalarSolutionStepVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteScalarSolutionStepVariable.md)
- [WriteVectorSolutionStepVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteVectorSolutionStepVariable.md)
- [WriteFlagContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable.md)
- [WriteScalarContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable.md)
- [WriteIntegrationScalarContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable.md)
- [WriteVectorContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable.md)
- [WriteIntegrationVectorContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable.md)
- [WritePropertiesIdsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WritePropertiesIdsToFile.md)
- [WriteIdsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteIdsToFile.md)