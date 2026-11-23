---
type: C++ Method
title: WriteIntegrationScalarContainerVariable
resource: input_output/vtk_output.cpp#L920-L954
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
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
  - target: functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerIntegrationResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteIntegrationScalarContainerVariable( const TContainerType& rContainer, const Variable<TVarType>& rVariable, std::ofstream& rFileStream) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [WriteScalarDataToFile](../../../../functions/input_output/vtk_output/WriteScalarDataToFile.md)

# Called by

- [WriteGeometricalContainerIntegrationResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerIntegrationResults.md)