---
type: C++ Method
title: Allocate
resource: includes/thread_fixed_size_memory_pool.h#L101-L123
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/AddChunk
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/chunk/Chunk/IsFull
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/chunk/Chunk/IsReleased
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void* Allocate()`

# Calls

- [AddChunk](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/AddChunk.md)
- [IsFull](../../../../functions/includes/chunk/Chunk/IsFull.md)
- [IsReleased](../../../../functions/includes/chunk/Chunk/IsReleased.md)