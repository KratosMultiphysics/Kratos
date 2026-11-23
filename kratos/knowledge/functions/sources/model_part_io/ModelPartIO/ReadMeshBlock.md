---
type: C++ Method
title: ReadMeshBlock
resource: sources/model_part_io.cpp#L4109-L4179
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/NumberOfMeshes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadMeshBlock(ModelPart& rModelPart)`

# Calls

- [NumberOfMeshes](../../../../functions/includes/model_part/NumberOfMeshes.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReadMeshDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshDataBlock.md)
- [ReadMeshNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshNodesBlock.md)
- [ReadMeshElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshElementsBlock.md)
- [ReadMeshConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshConditionsBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)

# Called by

- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)