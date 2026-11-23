---
type: C++ Method
title: ReadConstraintVectorialVariableData
resource: sources/model_part_io.cpp#L3323-L3352
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadConstraintVectorialVariableData(MasterSlaveConstraintContainerType& rThisConstraints, const TVariableType& rVariable, TDataType Dummy)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedConstraintId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId.md)

# Called by

- [ReadConstraintDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintDataBlock.md)