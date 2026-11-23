---
type: C++ Function
title: filename
resource: sources/kratos_filesystem.cpp#L79-L82
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string filename(const std::string& rPath)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE.md)
- [FileNameDataCollector](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector.md)
- [GetFileNameDataList](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList.md)