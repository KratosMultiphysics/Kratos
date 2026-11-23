---
type: C++ Method
title: RetrieveFileNameData
resource: utilities/file_name_data_collector.cpp#L138-L162
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/file_name_data_collector/UpdateFileNameData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameData/SetFileName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_file_name_data_collector/check_file_names
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_SortListOfFileNameData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_file_name_data_collector/check_sorted_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool FileNameDataCollector::RetrieveFileNameData( FileNameData& rFileNameData, const std::string& rFileNameWithoutPath) const`

# Calls

- [UpdateFileNameData](../../../../functions/utilities/file_name_data_collector/UpdateFileNameData.md)
- [SetFileName](../../../../functions/utilities/file_name_data_collector/FileNameData/SetFileName.md)

# Called by

- [check_file_names](../../../../functions/tests/test_file_name_data_collector/check_file_names.md)
- [test_SortListOfFileNameData](../../../../functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_SortListOfFileNameData.md)
- [check_sorted_list](../../../../functions/tests/test_file_name_data_collector/check_sorted_list.md)
- [GetFileNameDataList](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList.md)