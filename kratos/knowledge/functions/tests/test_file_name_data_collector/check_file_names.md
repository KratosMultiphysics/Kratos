---
type: Python Function
title: check_file_names
resource: tests/test_file_name_data_collector.py#L68-L78
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/RetrieveFileNameData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameData/SetFileName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_RetrieveFileNameData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def check_file_names(list_of_matching_indices, pattern):`

# Calls

- [FileNameDataCollector](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector.md)
- [RetrieveFileNameData](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/RetrieveFileNameData.md)
- [SetFileName](../../../functions/utilities/file_name_data_collector/FileNameData/SetFileName.md)

# Called by

- [test_RetrieveFileNameData](../../../functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_RetrieveFileNameData.md)