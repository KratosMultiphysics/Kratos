---
type: C++ Method
title: ReadNodesNumber
resource: sources/model_part_io.cpp#L114-L134
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ResetInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CountNodesInBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t ModelPartIO::ReadNodesNumber()`

# Calls

- [ResetInput](../../../../functions/sources/model_part_io/ModelPartIO/ResetInput.md)
- [CountNodesInBlock](../../../../functions/sources/model_part_io/ModelPartIO/CountNodesInBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)