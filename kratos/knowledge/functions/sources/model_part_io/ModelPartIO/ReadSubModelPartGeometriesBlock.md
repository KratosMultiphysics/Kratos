---
type: C++ Method
title: ReadSubModelPartGeometriesBlock
resource: sources/model_part_io.cpp#L4671-L4694
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadSubModelPartGeometriesBlock( ModelPart& rMainModelPart, ModelPart& rSubModelPart)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock.md)