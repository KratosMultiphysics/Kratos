---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/sources/test_filesystem.cpp#L29-L74
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/parent_path
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/filename
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/parallel_utilities/IndexPartition/IndexPartition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/ResolveSymlinks
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(ListDirectory, KratosCoreFastSuite)`

# Calls

- [parent_path](../../../../../functions/sources/kratos_filesystem/parent_path.md)
- [filename](../../../../../functions/sources/kratos_filesystem/filename.md)
- [IndexPartition](../../../../../functions/utilities/parallel_utilities/IndexPartition/IndexPartition.md)
- [MPISafeCreateDirectories](../../../../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)
- [emplace_back](../../../../../functions/containers/pointer_vector/emplace_back.md)
- [ResolveSymlinks](../../../../../functions/sources/kratos_filesystem/FilesystemExtensions/ResolveSymlinks.md)