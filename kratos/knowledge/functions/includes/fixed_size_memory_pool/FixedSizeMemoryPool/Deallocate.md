---
type: C++ Method
title: Deallocate
resource: includes/fixed_size_memory_pool.h#L89-L103
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/openmp_utils/GetCurrentNumberOfThreads
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Deallocate(void* pPointrerToRelease)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)
- [GetCurrentNumberOfThreads](../../../../functions/utilities/openmp_utils/GetCurrentNumberOfThreads.md)