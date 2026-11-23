---
type: C++ Method
title: ReadConditionalDataBlock
resource: sources/model_part_io.cpp#L3065-L3095
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionalScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionalVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadInitialValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadConditionalDataBlock(ConditionsContainerType& rThisConditions)`

# Calls

- [ReadConditionalScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionalScalarVariableData.md)
- [ReadConditionalVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionalVectorialVariableData.md)

# Called by

- [ReadInitialValues](../../../../functions/sources/model_part_io/ModelPartIO/ReadInitialValues.md)
- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)