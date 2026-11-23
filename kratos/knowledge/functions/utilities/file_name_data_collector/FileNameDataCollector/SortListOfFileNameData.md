---
type: C++ Method
title: SortListOfFileNameData
resource: utilities/file_name_data_collector.cpp#L204-L247
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/file_name_data_collector/FileNameData/GetStep
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_other_utilities_to_python/GetSortedListOfFileNameDataHelper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetSortedFileNamesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FileNameDataCollector::SortListOfFileNameData( std::vector<FileNameData>& rFileNameDataList, const std::vector<std::string>& rFlagsSortingOrder)`

# Calls

- [GetStep](../../../../functions/utilities/file_name_data_collector/FileNameData/GetStep.md)

# Called by

- [GetSortedListOfFileNameDataHelper](../../../../functions/python/add_other_utilities_to_python/GetSortedListOfFileNameDataHelper.md)
- [GetSortedFileNamesList](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetSortedFileNamesList.md)