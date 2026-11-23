---
type: C++ Method
title: WriteIntegrationVectorContainerVariable
resource: input_output/vtk_output.cpp#L986-L1028
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/WriteVectorDataToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerIntegrationResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteIntegrationVectorContainerVariable( const TContainerType& rContainer, const Variable<TVarType>& rVariable, std::ofstream& rFileStream) const`

# Calls

- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [WriteVectorDataToFile](../../../../functions/input_output/vtk_output/WriteVectorDataToFile.md)

# Called by

- [WriteGeometricalContainerIntegrationResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerIntegrationResults.md)