---
type: C++ Method
title: GetNumberOfAllocatedChunks
resource: includes/fixed_size_memory_pool.h#L119-L124
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/GetCurrentNumberOfThreads
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/GetNumberOfChunks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/GetNumberOfReleasedChunks
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_memory_pool/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t GetNumberOfAllocatedChunks() const`

# Calls

- [GetCurrentNumberOfThreads](../../../../functions/utilities/openmp_utils/GetCurrentNumberOfThreads.md)
- [GetNumberOfChunks](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/GetNumberOfChunks.md)
- [GetNumberOfReleasedChunks](../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/GetNumberOfReleasedChunks.md)

# Called by

- [PrintData](../../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/PrintData.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_memory_pool/KRATOS_TEST_CASE_IN_SUITE.md)