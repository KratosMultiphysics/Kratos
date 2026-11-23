---
type: C++ Method
title: GetNumberOfBlocks
resource: includes/chunk.h#L191-L193
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/chunk/Chunk/AllocatableDataSize
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/chunk/Chunk/MemoryOverhead
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/chunk/Chunk/IsEmpty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_chunk/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType GetNumberOfBlocks() const`

# Calls

- [AllocatableDataSize](../../../../functions/includes/chunk/Chunk/AllocatableDataSize.md)

# Called by

- [MemoryOverhead](../../../../functions/includes/chunk/Chunk/MemoryOverhead.md)
- [IsEmpty](../../../../functions/includes/chunk/Chunk/IsEmpty.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_chunk/KRATOS_TEST_CASE_IN_SUITE.md)