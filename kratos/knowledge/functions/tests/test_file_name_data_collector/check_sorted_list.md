---
type: Python Function
title: check_sorted_list
resource: tests/test_file_name_data_collector.py#L114-L118
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/RetrieveFileNameData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetSortedFileNamesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def check_sorted_list(sorting_order, sorted_indices_list):`

# Calls

- [RetrieveFileNameData](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/RetrieveFileNameData.md)
- [GetSortedFileNamesList](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetSortedFileNamesList.md)