---
type: Python Method
title: test_GetSortedFileNamesList
resource: tests/test_file_name_data_collector.py#L124-L161
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_GetSortedFileNamesList(self):`

# Calls

- [DeleteDirectoryIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)
- [FileNameDataCollector](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector.md)