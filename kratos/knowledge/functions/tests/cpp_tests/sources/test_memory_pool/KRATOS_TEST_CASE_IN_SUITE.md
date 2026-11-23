---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/sources/test_memory_pool.cpp#L29-L58
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/GetNumberOfChunks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/GetNumberOfAllocatedChunks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/openmp_utils/GetCurrentNumberOfThreads
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/exception/Exception/what
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Timer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(ThreadFixedSizeMemoryPoolAllocationDeallocation, KratosCoreFastSuite)`

# Calls

- [GetNumberOfChunks](../../../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/GetNumberOfChunks.md)
- [GetNumberOfAllocatedChunks](../../../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/GetNumberOfAllocatedChunks.md)
- [GetCurrentNumberOfThreads](../../../../../functions/utilities/openmp_utils/GetCurrentNumberOfThreads.md)
- [what](../../../../../functions/sources/exception/Exception/what.md)
- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)
- [Timer](../../../../../functions/utilities/timer/Timer/Timer.md)