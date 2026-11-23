---
type: C++ Method
title: Deallocate
resource: includes/thread_fixed_size_memory_pool.h#L125-L151
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/DeallocateFromAvailableChunk
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/chunk/Chunk/IsFull
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/DeallocateFromFullChunk
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Deallocate(void* pPointerToRelease)`

# Calls

- [DeallocateFromAvailableChunk](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/DeallocateFromAvailableChunk.md)
- [IsFull](../../../../functions/includes/chunk/Chunk/IsFull.md)
- [DeallocateFromFullChunk](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/DeallocateFromFullChunk.md)