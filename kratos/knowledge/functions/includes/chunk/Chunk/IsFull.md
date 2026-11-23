---
type: C++ Method
title: IsFull
resource: includes/chunk.h#L230-L232
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/chunk/Chunk/Allocate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/Allocate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/Deallocate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsFull()`

# Called by

- [Allocate](../../../../functions/includes/chunk/Chunk/Allocate.md)
- [Allocate](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/Allocate.md)
- [Deallocate](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/Deallocate.md)