---
type: C++ Method
title: GetNumberOfChunks
resource: includes/thread_fixed_size_memory_pool.h#L191-L193
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/GetNumberOfAllocatedChunks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_memory_pool/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t GetNumberOfChunks() const`

# Called by

- [GetNumberOfAllocatedChunks](../../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/GetNumberOfAllocatedChunks.md)
- [PrintData](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/PrintData.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_memory_pool/KRATOS_TEST_CASE_IN_SUITE.md)