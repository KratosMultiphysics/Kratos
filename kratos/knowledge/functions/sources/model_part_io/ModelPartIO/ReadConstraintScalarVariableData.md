---
type: C++ Method
title: ReadConstraintScalarVariableData
resource: sources/model_part_io.cpp#L3290-L3320
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

`void ModelPartIO::ReadConstraintScalarVariableData(MasterSlaveConstraintContainerType& rThisConstraints, const TVariableType& rVariable)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedConstraintId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId.md)

# Called by

- [ReadConstraintDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintDataBlock.md)