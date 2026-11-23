---
type: C++ Method
title: CheckEndBlock
resource: sources/model_part_io.cpp#L1591-L1601
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckStatement
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPartDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadTableBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CountNodesInBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadPropertiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalDofVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementalScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementalVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionalScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionalVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometryScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometryVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshPropertiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartPropertiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideFlagVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideDofVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ScanNodeBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartElementsAndConditionsIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartEntitiesIds
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ModelPartIO::CheckEndBlock(std::string const& BlockName, std::string& rWord)`

# Calls

- [CheckStatement](../../../../functions/sources/model_part_io/ModelPartIO/CheckStatement.md)

# Called by

- [FillNodalConnectivitiesFromGeometryBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlockInList.md)
- [FillNodalConnectivitiesFromElementBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlockInList.md)
- [FillNodalConnectivitiesFromConditionBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlockInList.md)
- [FillNodalConnectivitiesFromConstraintBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList.md)
- [ReadModelPartDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPartDataBlock.md)
- [ReadTableBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadTableBlock.md)
- [ReadNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodesBlock.md)
- [CountNodesInBlock](../../../../functions/sources/model_part_io/ModelPartIO/CountNodesInBlock.md)
- [ReadPropertiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadPropertiesBlock.md)
- [ReadGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesBlock.md)
- [ReadElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsBlock.md)
- [ReadConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionsBlock.md)
- [ReadConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock.md)
- [ReadNodalDofVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDofVariableData.md)
- [ReadNodalFlags](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalFlags.md)
- [ReadNodalScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalScalarVariableData.md)
- [ReadNodalVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalVectorialVariableData.md)
- [ReadElementalScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementalScalarVariableData.md)
- [ReadElementalVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementalVectorialVariableData.md)
- [ReadConditionalScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionalScalarVariableData.md)
- [ReadConditionalVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionalVectorialVariableData.md)
- [ReadGeometryScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometryScalarVariableData.md)
- [ReadGeometryVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometryVectorialVariableData.md)
- [ReadConstraintScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintScalarVariableData.md)
- [ReadConstraintVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintVectorialVariableData.md)
- [ReadGeometriesConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivitiesBlock.md)
- [ReadElementsConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivitiesBlock.md)
- [ReadConditionsConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivitiesBlock.md)
- [ReadConstraintsConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock.md)
- [FillNodalConnectivitiesFromGeometryBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlock.md)
- [FillNodalConnectivitiesFromElementBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlock.md)
- [FillNodalConnectivitiesFromConditionBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlock.md)
- [FillNodalConnectivitiesFromConstraintBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlock.md)
- [ReadCommunicatorDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock.md)
- [ReadCommunicatorLocalNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock.md)
- [ReadCommunicatorGhostNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock.md)
- [ReadMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshBlock.md)
- [ReadMeshDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshDataBlock.md)
- [ReadMeshNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshNodesBlock.md)
- [ReadMeshElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshElementsBlock.md)
- [ReadMeshConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshConditionsBlock.md)
- [ReadMeshPropertiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshPropertiesBlock.md)
- [ReadSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock.md)
- [ReadSubModelPartTablesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock.md)
- [ReadSubModelPartPropertiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartPropertiesBlock.md)
- [ReadSubModelPartNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartNodesBlock.md)
- [ReadSubModelPartElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartElementsBlock.md)
- [ReadSubModelPartConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConditionsBlock.md)
- [ReadSubModelPartGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartGeometriesBlock.md)
- [ReadSubModelPartConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConstraintsBlock.md)
- [DivideNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideNodesBlock.md)
- [DivideGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideGeometriesBlock.md)
- [DivideElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideElementsBlock.md)
- [DivideConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConditionsBlock.md)
- [DivideConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock.md)
- [DivideFlagVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideFlagVariableData.md)
- [DivideDofVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideDofVariableData.md)
- [DivideVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideVectorialVariableData.md)
- [DivideScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData.md)
- [DivideMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshBlock.md)
- [DivideSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock.md)
- [DivideMeshNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshNodesBlock.md)
- [DivideMeshElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshElementsBlock.md)
- [DivideMeshConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshConditionsBlock.md)
- [DivideMeshConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshConstraintsBlock.md)
- [DivideSubModelPartNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartNodesBlock.md)
- [DivideSubModelPartElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartElementsBlock.md)
- [DivideSubModelPartConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConditionsBlock.md)
- [DivideSubModelPartConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConstraintsBlock.md)
- [DivideSubModelPartGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartGeometriesBlock.md)
- [ScanNodeBlock](../../../../functions/sources/model_part_io/ModelPartIO/ScanNodeBlock.md)
- [ReadSubModelPartElementsAndConditionsIds](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartElementsAndConditionsIds.md)
- [ReadSubModelPartEntitiesIds](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartEntitiesIds.md)