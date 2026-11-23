---
type: C++ Method
title: ResolveSymlinks
resource: sources/kratos_filesystem.cpp#L135-L151
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/kratos_parameters/Parameters/SolveIncludes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::filesystem::path FilesystemExtensions::ResolveSymlinks(const std::filesystem::path& rPath)`

# Called by

- [SolveIncludes](../../../../functions/sources/kratos_parameters/Parameters/SolveIncludes.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE.md)