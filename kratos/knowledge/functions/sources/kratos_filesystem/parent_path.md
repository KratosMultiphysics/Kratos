---
type: C++ Function
title: parent_path
resource: sources/kratos_filesystem.cpp#L74-L77
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtu_output/GetFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string parent_path(const std::string& rPath)`

# Called by

- [GetFileName](../../../functions/input_output/vtu_output/GetFileName.md)
- [PrintOutput](../../../functions/input_output/vtu_output/VtuOutput/PrintOutput.md)
- [DivideInputToPartitions](../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitions.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE.md)
- [FileNameDataCollector](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector.md)