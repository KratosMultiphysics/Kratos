---
type: C++ Method
title: AllocatableDataSize
resource: includes/chunk.h#L299-L305
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/chunk/Chunk/GetHeaderSize
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/chunk/Chunk/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/chunk/Chunk/GetNumberOfBlocks
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType AllocatableDataSize() const`

# Calls

- [GetHeaderSize](../../../../functions/includes/chunk/Chunk/GetHeaderSize.md)

# Called by

- [Initialize](../../../../functions/includes/chunk/Chunk/Initialize.md)
- [GetNumberOfBlocks](../../../../functions/includes/chunk/Chunk/GetNumberOfBlocks.md)