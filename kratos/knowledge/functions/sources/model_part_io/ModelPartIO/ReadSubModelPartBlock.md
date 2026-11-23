---
type: C++ Method
title: ReadSubModelPartBlock
resource: sources/model_part_io.cpp#L4426-L4483
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartPropertiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadSubModelPartBlock(ModelPart& rMainModelPart, ModelPart& rParentModelPart)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [ReadSubModelPartDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartDataBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)
- [ReadSubModelPartTablesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock.md)
- [ReadSubModelPartPropertiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartPropertiesBlock.md)
- [ReadSubModelPartNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartNodesBlock.md)
- [ReadSubModelPartElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartElementsBlock.md)
- [ReadSubModelPartConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConditionsBlock.md)
- [ReadSubModelPartGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartGeometriesBlock.md)
- [ReadSubModelPartConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartConstraintsBlock.md)

# Called by

- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)