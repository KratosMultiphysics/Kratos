---
type: C++ Method
title: Deallocate
resource: includes/memory_pool.h#L81-L83
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/memory_pool/MemoryPool/GetPoolWithBlockSize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void Deallocate(void* pPointrerToRelease, std::size_t ObjectSizeInBytes)`

# Calls

- [GetPoolWithBlockSize](../../../../functions/includes/memory_pool/MemoryPool/GetPoolWithBlockSize.md)