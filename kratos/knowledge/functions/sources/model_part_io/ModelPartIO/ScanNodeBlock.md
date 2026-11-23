---
type: C++ Method
title: ScanNodeBlock
resource: sources/model_part_io.cpp#L6579-L6597
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ScanNodeBlock()`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadNodalGraph](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraph.md)
- [ReadNodalGraphFromEntitiesList](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList.md)