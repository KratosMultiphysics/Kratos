---
type: C++ Function
title: WriteVectorDataToFile
resource: input_output/vtk_output.h#L462-L475
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/ForceBigEndian
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodesToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteVectorSolutionStepVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteVectorDataToFile(const TData& rData, std::ofstream& rFileStream) const`

# Calls

- [ForceBigEndian](../../../functions/input_output/vtk_output/VtkOutput/ForceBigEndian.md)

# Called by

- [WriteNodesToFile](../../../functions/input_output/vtk_output/VtkOutput/WriteNodesToFile.md)
- [WriteVectorSolutionStepVariable](../../../functions/input_output/vtk_output/VtkOutput/WriteVectorSolutionStepVariable.md)
- [WriteVectorContainerVariable](../../../functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable.md)
- [WriteIntegrationVectorContainerVariable](../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable.md)