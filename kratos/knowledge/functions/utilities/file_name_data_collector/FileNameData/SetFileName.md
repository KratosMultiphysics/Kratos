---
type: C++ Method
title: SetFileName
resource: utilities/file_name_data_collector.h#L80
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_file_name_data_collector/check_file_names
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/RetrieveFileNameData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetFileName(const std::string& rFileName)`

# Called by

- [check_file_names](../../../../functions/tests/test_file_name_data_collector/check_file_names.md)
- [RetrieveFileNameData](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/RetrieveFileNameData.md)