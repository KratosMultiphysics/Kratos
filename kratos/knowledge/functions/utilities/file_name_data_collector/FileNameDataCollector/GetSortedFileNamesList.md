---
type: C++ Method
title: GetSortedFileNamesList
resource: utilities/file_name_data_collector.cpp#L185-L202
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/SortListOfFileNameData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_file_name_data_collector/check_sorted_list
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> FileNameDataCollector::GetSortedFileNamesList( const std::vector<std::string>& rFlagsSortingOrder) const`

# Calls

- [GetFileNameDataList](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList.md)
- [SortListOfFileNameData](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/SortListOfFileNameData.md)

# Called by

- [check_sorted_list](../../../../functions/tests/test_file_name_data_collector/check_sorted_list.md)