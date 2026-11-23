---
type: C++ Method
title: SkippableEntity
resource: utilities/input_output_utilities.cpp#L29-L82
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/CollectPartData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkCellListSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteConnectivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteCellType
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

`bool InputOutputUtilities::SkippableEntity( const TEntityType& rEntity, const std::string& rWarningLabel )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CollectPartData](../../../../functions/input_output/ensight_output/EnSightOutput/CollectPartData.md)
- [DetermineVtkContainerSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [DetermineVtkCellListSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkCellListSize.md)
- [WriteConnectivity](../../../../functions/input_output/vtk_output/VtkOutput/WriteConnectivity.md)
- [WriteCellType](../../../../functions/input_output/vtk_output/VtkOutput/WriteCellType.md)
- [WriteFlagContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable.md)
- [WriteScalarContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable.md)
- [WriteIntegrationScalarContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable.md)
- [WriteVectorContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable.md)
- [WriteIntegrationVectorContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable.md)
- [WritePropertiesIdsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WritePropertiesIdsToFile.md)
- [WriteIdsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteIdsToFile.md)