---
type: C++ Method
title: ReadModelPartDataBlock
resource: sources/model_part_io.cpp#L1612-L1692
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadModelPartDataBlock(ModelPart& rModelPart, const bool IsSubmodelpart)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)