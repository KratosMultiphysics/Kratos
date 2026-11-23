---
type: C++ Method
title: ReadConstraintDataBlock
resource: sources/model_part_io.cpp#L3257-L3287
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintVectorialVariableData
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

`void ModelPartIO::ReadConstraintDataBlock(MasterSlaveConstraintContainerType& rThisConstraints)`

# Calls

- [ReadConstraintScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintScalarVariableData.md)
- [ReadConstraintVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintVectorialVariableData.md)

# Called by

- [ReadInitialValues](../../../../functions/sources/model_part_io/ModelPartIO/ReadInitialValues.md)
- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)