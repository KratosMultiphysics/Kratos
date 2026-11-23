---
type: C++ Function
title: GetCurrentNumberOfThreads
resource: utilities/openmp_utils.h#L73-L80
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/FixedSizeMemoryPool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/Deallocate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/MemoryUsed
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/MemoryOverhead
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/GetNumberOfAllocatedChunks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_memory_pool/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static int GetCurrentNumberOfThreads()`

# Called by

- [FixedSizeMemoryPool](../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/FixedSizeMemoryPool.md)
- [Deallocate](../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/Deallocate.md)
- [MemoryUsed](../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/MemoryUsed.md)
- [MemoryOverhead](../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/MemoryOverhead.md)
- [GetNumberOfAllocatedChunks](../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/GetNumberOfAllocatedChunks.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_memory_pool/KRATOS_TEST_CASE_IN_SUITE.md)