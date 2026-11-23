---
type: C++ Method
title: DeallocateFromFullChunk
resource: includes/thread_fixed_size_memory_pool.h#L281-L288
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/ReleaseChunk
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/Deallocate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DeallocateFromFullChunk(void* pPointrerToRelease, Chunk* pChunk)`

# Calls

- [ReleaseChunk](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/ReleaseChunk.md)

# Called by

- [Deallocate](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/Deallocate.md)