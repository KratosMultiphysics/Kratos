---
type: C++ Method
title: ReadConditionalScalarVariableData
resource: sources/model_part_io.cpp#L3098-L3127
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionalDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadConditionalScalarVariableData(ConditionsContainerType& rThisConditions, const TVariableType& rVariable)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadConditionalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionalDataBlock.md)