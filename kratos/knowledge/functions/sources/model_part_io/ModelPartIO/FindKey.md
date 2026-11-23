---
type: C++ Method
title: FindKey
resource: sources/model_part_io.cpp#L6356-L6368
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
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
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock
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
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartPropertiesBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename TContainerType::iterator ModelPartIO::FindKey(TContainerType& ThisContainer , TKeyType ThisKey, std::string ComponentName)`

# Called by

- [ReadGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesBlock.md)
- [ReadElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsBlock.md)
- [ReadConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionsBlock.md)
- [ReadConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock.md)
- [ReadNodalDofVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDofVariableData.md)
- [ReadNodalFlags](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalFlags.md)
- [ReadNodalScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalScalarVariableData.md)
- [ReadNodalVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalVectorialVariableData.md)
- [ReadCommunicatorLocalNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock.md)
- [ReadCommunicatorGhostNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock.md)
- [ReadMeshNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshNodesBlock.md)
- [ReadMeshElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshElementsBlock.md)
- [ReadMeshConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshConditionsBlock.md)
- [ReadSubModelPartTablesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock.md)
- [ReadSubModelPartPropertiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartPropertiesBlock.md)