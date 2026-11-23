---
type: C++ Method
title: CountNodesInBlock
resource: sources/model_part_io.cpp#L2023-L2064
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodesNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t ModelPartIO::CountNodesInBlock()`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadNodesNumber](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodesNumber.md)