---
type: C++ Function
title: WriteScalarDataToFile
resource: input_output/vtk_output.h#L444-L453
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/ForceBigEndian
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteConnectivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteCellType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteScalarSolutionStepVariable
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
  - target: functions/input_output/vtk_output/VtkOutput/WritePropertiesIdsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIdsToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteScalarDataToFile(const TData& rData, std::ofstream& rFileStream) const`

# Calls

- [ForceBigEndian](../../../functions/input_output/vtk_output/VtkOutput/ForceBigEndian.md)

# Called by

- [WriteConnectivity](../../../functions/input_output/vtk_output/VtkOutput/WriteConnectivity.md)
- [WriteCellType](../../../functions/input_output/vtk_output/VtkOutput/WriteCellType.md)
- [WriteScalarSolutionStepVariable](../../../functions/input_output/vtk_output/VtkOutput/WriteScalarSolutionStepVariable.md)
- [WriteFlagContainerVariable](../../../functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable.md)
- [WriteScalarContainerVariable](../../../functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable.md)
- [WriteIntegrationScalarContainerVariable](../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable.md)
- [WritePropertiesIdsToFile](../../../functions/input_output/vtk_output/VtkOutput/WritePropertiesIdsToFile.md)
- [WriteIdsToFile](../../../functions/input_output/vtk_output/VtkOutput/WriteIdsToFile.md)