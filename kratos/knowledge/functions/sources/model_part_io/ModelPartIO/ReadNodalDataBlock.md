---
type: C++ Method
title: ReadNodalDataBlock
resource: sources/model_part_io.cpp#L2590-L2656
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalDofVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalVectorialVariableData
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

`void ModelPartIO::ReadNodalDataBlock(ModelPart& rThisModelPart)`

# Calls

- [ReadNodalFlags](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalFlags.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)
- [ReadNodalScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalScalarVariableData.md)
- [ReadNodalDofVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDofVariableData.md)
- [ReadNodalVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalVectorialVariableData.md)

# Called by

- [ReadInitialValues](../../../../functions/sources/model_part_io/ModelPartIO/ReadInitialValues.md)
- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)