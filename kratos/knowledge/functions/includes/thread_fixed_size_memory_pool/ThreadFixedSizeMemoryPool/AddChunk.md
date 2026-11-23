---
type: C++ Method
title: AddChunk
resource: includes/thread_fixed_size_memory_pool.h#L262-L272
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/Allocate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddChunk()`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)
- [emplace_back](../../../../functions/containers/pointer_vector/emplace_back.md)

# Called by

- [Allocate](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/Allocate.md)