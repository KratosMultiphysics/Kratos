---
type: C++ Module
title: vtk_output
resource: input_output/vtk_output.cpp#L1-L1181
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/map
    resolved_by: tree-sitter
    confidence: exact
  - target: external/vtk-output-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-model-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-filesystem-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-fast-transfer-between-model-parts-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-input-output-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/input-output-vtk-definitions-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetDefaultParameters](../../functions/input_output/vtk_output/VtkOutput/GetDefaultParameters.md)
- [VtkOutput](../../functions/input_output/vtk_output/VtkOutput/VtkOutput.md)
- [GetEntityType](../../functions/input_output/vtk_output/VtkOutput/GetEntityType.md)
- [PrepareGaussPointResults](../../functions/input_output/vtk_output/VtkOutput/PrepareGaussPointResults.md)
- [PrintOutput](../../functions/input_output/vtk_output/VtkOutput/PrintOutput.md)
- [WriteModelPartToFile](../../functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile.md)
- [GetOutputFileName](../../functions/input_output/vtk_output/VtkOutput/GetOutputFileName.md)
- [Initialize](../../functions/input_output/vtk_output/VtkOutput/Initialize.md)
- [CreateMapFromKratosIdToVTKId](../../functions/input_output/vtk_output/VtkOutput/CreateMapFromKratosIdToVTKId.md)
- [WriteHeaderToFile](../../functions/input_output/vtk_output/VtkOutput/WriteHeaderToFile.md)
- [WriteMeshToFile](../../functions/input_output/vtk_output/VtkOutput/WriteMeshToFile.md)
- [WriteNodesToFile](../../functions/input_output/vtk_output/VtkOutput/WriteNodesToFile.md)
- [WriteConditionsAndElementsToFile](../../functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile.md)
- [DetermineVtkContainerSize](../../functions/input_output/vtk_output/VtkOutput/DetermineVtkContainerSize.md)
- [DetermineVtkCellListSize](../../functions/input_output/vtk_output/VtkOutput/DetermineVtkCellListSize.md)
- [WriteConnectivity](../../functions/input_output/vtk_output/VtkOutput/WriteConnectivity.md)
- [WriteCellType](../../functions/input_output/vtk_output/VtkOutput/WriteCellType.md)
- [IsCompatibleVariable](../../functions/input_output/vtk_output/VtkOutput/IsCompatibleVariable.md)
- [WriteNodalResultsToFile](../../functions/input_output/vtk_output/VtkOutput/WriteNodalResultsToFile.md)
- [WriteElementResultsToFile](../../functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile.md)
- [WriteConditionResultsToFile](../../functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile.md)
- [WriteNodalContainerResults](../../functions/input_output/vtk_output/VtkOutput/WriteNodalContainerResults.md)
- [WriteGeometricalContainerResults](../../functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerResults.md)
- [WriteGeometricalContainerIntegrationResults](../../functions/input_output/vtk_output/VtkOutput/WriteGeometricalContainerIntegrationResults.md)
- [WriteNodalScalarValues](../../functions/input_output/vtk_output/VtkOutput/WriteNodalScalarValues.md)
- [WriteNodalVectorValues](../../functions/input_output/vtk_output/VtkOutput/WriteNodalVectorValues.md)
- [WriteScalarSolutionStepVariable](../../functions/input_output/vtk_output/VtkOutput/WriteScalarSolutionStepVariable.md)
- [WriteVectorSolutionStepVariable](../../functions/input_output/vtk_output/VtkOutput/WriteVectorSolutionStepVariable.md)
- [WriteFlagContainerVariable](../../functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable.md)
- [WriteScalarContainerVariable](../../functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable.md)
- [WriteIntegrationScalarContainerVariable](../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable.md)
- [WriteVectorContainerVariable](../../functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable.md)
- [WriteIntegrationVectorContainerVariable](../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable.md)
- [ForceBigEndian](../../functions/input_output/vtk_output/VtkOutput/ForceBigEndian.md)
- [ReorderConnectivity](../../functions/input_output/vtk_output/VtkOutput/ReorderConnectivity.md)
- [WritePropertiesIdsToFile](../../functions/input_output/vtk_output/VtkOutput/WritePropertiesIdsToFile.md)
- [WriteIdsToFile](../../functions/input_output/vtk_output/VtkOutput/WriteIdsToFile.md)
- [WriteModelPartWithoutNodesToFile](../../functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile.md)

# Imports

- `<map>`
- `"vtk_output.h"`
- `"containers/model.h"`
- `"includes/kratos_filesystem.h"`
- `"processes/fast_transfer_between_model_parts_process.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/reduction_utilities.h"`
- `"utilities/input_output_utilities.h"`
- `"input_output/vtk_definitions.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)