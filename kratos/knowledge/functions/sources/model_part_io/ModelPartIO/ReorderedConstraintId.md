---
type: C++ Method
title: ReorderedConstraintId
resource: sources/model_part_io.cpp#L6739-L6744
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartEntitiesIds
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPartIO::SizeType ModelPartIO::ReorderedConstraintId(ModelPartIO::SizeType ConstraintId)`

# Called by

- [FillNodalConnectivitiesFromConstraintBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList.md)
- [ReadConstraintScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintScalarVariableData.md)
- [ReadConstraintVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintVectorialVariableData.md)
- [ReadConstraintsConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock.md)
- [DivideConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock.md)
- [DivideVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideVectorialVariableData.md)
- [DivideScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData.md)
- [DivideMeshConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshConstraintsBlock.md)
- [DivideSubModelPartConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConstraintsBlock.md)
- [ReadSubModelPartEntitiesIds](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartEntitiesIds.md)