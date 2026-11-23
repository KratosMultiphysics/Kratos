---
type: C++ Method
title: ReadInitialValues
resource: sources/model_part_io.cpp#L784-L815
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ResetInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometryDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadInitialValues(ModelPart& rThisModelPart)`

# Calls

- [ResetInput](../../../../functions/sources/model_part_io/ModelPartIO/ResetInput.md)
- [ReadNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock.md)
- [ReadElementalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementalDataBlock.md)
- [ReadConditionalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionalDataBlock.md)
- [ReadConstraintDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintDataBlock.md)
- [ReadGeometryDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometryDataBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)