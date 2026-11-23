---
type: C++ Module
title: model_part_io
resource: sources/model_part_io.cpp#L1-L6745
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/set
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-io-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-filesystem-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/input-output-logger-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-compare-elements-and-conditions-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-openmp-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-quaternion-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-string-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-timer-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ModelPartIO](../../functions/sources/model_part_io/ModelPartIO/ModelPartIO.md)
- [ModelPartIO](../../functions/sources/model_part_io/ModelPartIO/ModelPartIO-2.md)
- [ReadNode](../../functions/sources/model_part_io/ModelPartIO/ReadNode.md)
- [ReadNodes](../../functions/sources/model_part_io/ModelPartIO/ReadNodes.md)
- [ReadNodesNumber](../../functions/sources/model_part_io/ModelPartIO/ReadNodesNumber.md)
- [WriteNodes](../../functions/sources/model_part_io/ModelPartIO/WriteNodes.md)
- [ReadProperties](../../functions/sources/model_part_io/ModelPartIO/ReadProperties.md)
- [ReadProperties](../../functions/sources/model_part_io/ModelPartIO/ReadProperties-2.md)
- [WriteProperties](../../functions/sources/model_part_io/ModelPartIO/WriteProperties.md)
- [ReadGeometry](../../functions/sources/model_part_io/ModelPartIO/ReadGeometry.md)
- [ReadGeometries](../../functions/sources/model_part_io/ModelPartIO/ReadGeometries.md)
- [ReadGeometriesConnectivities](../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivities.md)
- [WriteGeometries](../../functions/sources/model_part_io/ModelPartIO/WriteGeometries.md)
- [ReadElement](../../functions/sources/model_part_io/ModelPartIO/ReadElement.md)
- [ReadElements](../../functions/sources/model_part_io/ModelPartIO/ReadElements.md)
- [ReadElementsConnectivities](../../functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivities.md)
- [WriteElements](../../functions/sources/model_part_io/ModelPartIO/WriteElements.md)
- [ReadConditions](../../functions/sources/model_part_io/ModelPartIO/ReadConditions.md)
- [ReadConditionsConnectivities](../../functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivities.md)
- [WriteConditions](../../functions/sources/model_part_io/ModelPartIO/WriteConditions.md)
- [ReadConstraints](../../functions/sources/model_part_io/ModelPartIO/ReadConstraints.md)
- [ReadConstraintsConnectivities](../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivities.md)
- [WriteConstraints](../../functions/sources/model_part_io/ModelPartIO/WriteConstraints.md)
- [ReadInitialValues](../../functions/sources/model_part_io/ModelPartIO/ReadInitialValues.md)
- [ReadMesh](../../functions/sources/model_part_io/ModelPartIO/ReadMesh.md)
- [WriteMesh](../../functions/sources/model_part_io/ModelPartIO/WriteMesh.md)
- [ReadModelPart](../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)
- [WriteModelPart](../../functions/sources/model_part_io/ModelPartIO/WriteModelPart.md)
- [ReadNodalGraph](../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraph.md)
- [ReadNodalGraphFromEntitiesList](../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList.md)
- [ReadNodalGraphFromEntitiesList](../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList-2.md)
- [FillNodalConnectivitiesFromGeometryBlockInList](../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlockInList.md)
- [FillNodalConnectivitiesFromElementBlockInList](../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlockInList.md)
- [FillNodalConnectivitiesFromConditionBlockInList](../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlockInList.md)
- [FillNodalConnectivitiesFromConstraintBlockInList](../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList.md)
- [DivideInputToPartitions](../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitions.md)
- [DivideInputToPartitions](../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitions-2.md)
- [DivideInputToPartitionsImpl](../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)
- [SkipBlock](../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)
- [CheckEndBlock](../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReadSubModelPartDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartDataBlock.md)
- [ReadModelPartDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadModelPartDataBlock.md)
- [WriteModelPartDataBlock](../../functions/sources/model_part_io/ModelPartIO/WriteModelPartDataBlock.md)
- [ReadTableBlock](../../functions/sources/model_part_io/ModelPartIO/ReadTableBlock.md)
- [ReadTableBlock](../../functions/sources/model_part_io/ModelPartIO/ReadTableBlock-2.md)
- [WriteTableBlock](../../functions/sources/model_part_io/ModelPartIO/WriteTableBlock.md)
- [WriteTableBlock](../../functions/sources/model_part_io/ModelPartIO/WriteTableBlock-2.md)
- [ReadNodesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadNodesBlock.md)
- [ReadNodesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadNodesBlock-2.md)
- [CountNodesInBlock](../../functions/sources/model_part_io/ModelPartIO/CountNodesInBlock.md)
- [ReadPropertiesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadPropertiesBlock.md)
- [ReadGeometriesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesBlock.md)
- [ReadGeometriesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesBlock-2.md)
- [ReadElementsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadElementsBlock.md)
- [ReadElementsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadElementsBlock-2.md)
- [ReadConditionsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadConditionsBlock.md)
- [ReadConditionsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadConditionsBlock-2.md)
- [ReadConstraintsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock.md)
- [ReadConstraintsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock-2.md)
- [ReadNodalDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock.md)
- [WriteNodalDataBlock](../../functions/sources/model_part_io/ModelPartIO/WriteNodalDataBlock.md)
- [WriteDataBlock](../../functions/sources/model_part_io/ModelPartIO/WriteDataBlock.md)
- [WriteDataBlock](../../functions/sources/model_part_io/ModelPartIO/WriteDataBlock-2.md)
- [ReadNodalDofVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadNodalDofVariableData.md)
- [ReadNodalFlags](../../functions/sources/model_part_io/ModelPartIO/ReadNodalFlags.md)
- [ReadNodalScalarVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadNodalScalarVariableData.md)
- [ReadNodalVectorialVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadNodalVectorialVariableData.md)
- [ReadElementalDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadElementalDataBlock.md)
- [ReadElementalScalarVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadElementalScalarVariableData.md)
- [ReadElementalVectorialVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadElementalVectorialVariableData.md)
- [ReadConditionalDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadConditionalDataBlock.md)
- [ReadConditionalScalarVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadConditionalScalarVariableData.md)
- [ReadConditionalVectorialVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadConditionalVectorialVariableData.md)
- [ReadGeometryDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadGeometryDataBlock.md)
- [ReadGeometryScalarVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadGeometryScalarVariableData.md)
- [ReadGeometryVectorialVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadGeometryVectorialVariableData.md)
- [ReadConstraintDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadConstraintDataBlock.md)
- [ReadConstraintScalarVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadConstraintScalarVariableData.md)
- [ReadConstraintVectorialVariableData](../../functions/sources/model_part_io/ModelPartIO/ReadConstraintVectorialVariableData.md)
- [ReadGeometriesConnectivitiesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivitiesBlock.md)
- [ReadElementsConnectivitiesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivitiesBlock.md)
- [ReadConditionsConnectivitiesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivitiesBlock.md)
- [ReadConstraintsConnectivitiesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock.md)
- [FillNodalConnectivitiesFromGeometryBlock](../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlock.md)
- [FillNodalConnectivitiesFromElementBlock](../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlock.md)
- [FillNodalConnectivitiesFromConditionBlock](../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlock.md)
- [FillNodalConnectivitiesFromConstraintBlock](../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlock.md)
- [ReadCommunicatorDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock.md)
- [ReadCommunicatorLocalNodesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock.md)
- [ReadCommunicatorGhostNodesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock.md)
- [ReadMeshBlock](../../functions/sources/model_part_io/ModelPartIO/ReadMeshBlock.md)
- [ReadMeshDataBlock](../../functions/sources/model_part_io/ModelPartIO/ReadMeshDataBlock.md)
- [ReadMeshNodesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadMeshNodesBlock.md)
- [ReadMeshElementsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadMeshElementsBlock.md)
- [ReadMeshConditionsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadMeshConditionsBlock.md)
- [ReadMeshPropertiesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadMeshPropertiesBlock.md)
- [ReadSubModelPartBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock.md)
- [WriteSubModelPartBlock](../../functions/sources/model_part_io/ModelPartIO/WriteSubModelPartBlock.md)
- [ReadSubModelPartTablesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock.md)
- [ReadSubModelPartPropertiesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartPropertiesBlock.md)
- [ReadSubModelPartNodesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartNodesBlock.md)
- [ReadSubModelPartElementsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartElementsBlock.md)
- [ReadSubModelPartConditionsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConditionsBlock.md)
- [ReadSubModelPartGeometriesBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartGeometriesBlock.md)
- [ReadSubModelPartConstraintsBlock](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConstraintsBlock.md)
- [DivideModelPartDataBlock](../../functions/sources/model_part_io/ModelPartIO/DivideModelPartDataBlock.md)
- [DivideTableBlock](../../functions/sources/model_part_io/ModelPartIO/DivideTableBlock.md)
- [DividePropertiesBlock](../../functions/sources/model_part_io/ModelPartIO/DividePropertiesBlock.md)
- [DivideNodesBlock](../../functions/sources/model_part_io/ModelPartIO/DivideNodesBlock.md)
- [DivideGeometriesBlock](../../functions/sources/model_part_io/ModelPartIO/DivideGeometriesBlock.md)
- [DivideElementsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideElementsBlock.md)
- [DivideConditionsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideConditionsBlock.md)
- [DivideConstraintsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock.md)
- [DivideNodalDataBlock](../../functions/sources/model_part_io/ModelPartIO/DivideNodalDataBlock.md)
- [DivideFlagVariableData](../../functions/sources/model_part_io/ModelPartIO/DivideFlagVariableData.md)
- [DivideDofVariableData](../../functions/sources/model_part_io/ModelPartIO/DivideDofVariableData.md)
- [DivideVectorialVariableData](../../functions/sources/model_part_io/ModelPartIO/DivideVectorialVariableData.md)
- [DivideElementalDataBlock](../../functions/sources/model_part_io/ModelPartIO/DivideElementalDataBlock.md)
- [DivideScalarVariableData](../../functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData.md)
- [DivideConditionalDataBlock](../../functions/sources/model_part_io/ModelPartIO/DivideConditionalDataBlock.md)
- [DivideConstraintDataBlock](../../functions/sources/model_part_io/ModelPartIO/DivideConstraintDataBlock.md)
- [DivideMeshBlock](../../functions/sources/model_part_io/ModelPartIO/DivideMeshBlock.md)
- [DivideSubModelPartBlock](../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock.md)
- [DivideMeshDataBlock](../../functions/sources/model_part_io/ModelPartIO/DivideMeshDataBlock.md)
- [DivideMeshNodesBlock](../../functions/sources/model_part_io/ModelPartIO/DivideMeshNodesBlock.md)
- [DivideMeshElementsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideMeshElementsBlock.md)
- [DivideMeshConditionsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideMeshConditionsBlock.md)
- [DivideMeshConstraintsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideMeshConstraintsBlock.md)
- [DivideSubModelPartDataBlock](../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartDataBlock.md)
- [DivideSubModelPartTableBlock](../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartTableBlock.md)
- [DivideSubModelPartNodesBlock](../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartNodesBlock.md)
- [DivideSubModelPartElementsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartElementsBlock.md)
- [DivideSubModelPartConditionsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConditionsBlock.md)
- [DivideSubModelPartConstraintsBlock](../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConstraintsBlock.md)
- [DivideSubModelPartGeometriesBlock](../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartGeometriesBlock.md)
- [WritePartitionIndices](../../functions/sources/model_part_io/ModelPartIO/WritePartitionIndices.md)
- [WriteCommunicatorData](../../functions/sources/model_part_io/ModelPartIO/WriteCommunicatorData.md)
- [WriteCommunicatorLocalNodes](../../functions/sources/model_part_io/ModelPartIO/WriteCommunicatorLocalNodes.md)
- [WriteInAllFiles](../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)
- [FindKey](../../functions/sources/model_part_io/ModelPartIO/FindKey.md)
- [ReadConstitutiveLawValue](../../functions/sources/model_part_io/ModelPartIO/ReadConstitutiveLawValue.md)
- [SkipWhiteSpaces](../../functions/sources/model_part_io/ModelPartIO/SkipWhiteSpaces.md)
- [GetCharacter](../../functions/sources/model_part_io/ModelPartIO/GetCharacter.md)
- [CheckStatement](../../functions/sources/model_part_io/ModelPartIO/CheckStatement.md)
- [ResetInput](../../functions/sources/model_part_io/ModelPartIO/ResetInput.md)
- [SwapStreamSource](../../functions/sources/model_part_io/ModelPartIO/SwapStreamSource.md)
- [CreatePartition](../../functions/sources/model_part_io/ModelPartIO/CreatePartition.md)
- [ScanNodeBlock](../../functions/sources/model_part_io/ModelPartIO/ScanNodeBlock.md)
- [ReadSubModelPartElementsAndConditionsIds](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartElementsAndConditionsIds.md)
- [ReadSubModelPartEntitiesIds](../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartEntitiesIds.md)
- [ReorderedNodeId](../../functions/sources/model_part_io/ModelPartIO/ReorderedNodeId.md)
- [ReorderedGeometryId](../../functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId.md)
- [ReorderedElementId](../../functions/sources/model_part_io/ModelPartIO/ReorderedElementId.md)
- [ReorderedConditionId](../../functions/sources/model_part_io/ModelPartIO/ReorderedConditionId.md)
- [ReorderedConstraintId](../../functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId.md)

# Imports

- `<set>`
- `"includes/model_part_io.h"`
- `"includes/kratos_filesystem.h"`
- `"input_output/logger.h"`
- `"utilities/compare_elements_and_conditions_utility.h"`
- `"utilities/openmp_utils.h"`
- `"utilities/quaternion.h"`
- `"utilities/string_utilities.h"`
- `"utilities/timer.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)