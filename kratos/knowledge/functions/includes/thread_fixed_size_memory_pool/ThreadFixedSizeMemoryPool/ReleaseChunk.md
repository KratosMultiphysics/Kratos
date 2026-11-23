---
type: C++ Method
title: ReleaseChunk
resource: includes/thread_fixed_size_memory_pool.h#L290-L294
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/DeallocateFromAvailableChunk
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/DeallocateFromFullChunk
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReleaseChunk(Chunk* pChunk)`

# Called by

- [DeallocateFromAvailableChunk](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/DeallocateFromAvailableChunk.md)
- [DeallocateFromFullChunk](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/DeallocateFromFullChunk.md)