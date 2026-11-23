---
type: C++ Method
title: GetPoolWithBlockSize
resource: includes/memory_pool.h#L99-L119
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/memory_pool/MemoryPool/Allocate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/memory_pool/MemoryPool/Deallocate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static FixedSizeMemoryPool* GetPoolWithBlockSize(std::size_t BlockSize)`

# Called by

- [Allocate](../../../../functions/includes/memory_pool/MemoryPool/Allocate.md)
- [Deallocate](../../../../functions/includes/memory_pool/MemoryPool/Deallocate.md)