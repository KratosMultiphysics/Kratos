---
type: C++ Function
title: GetSortedListOfFileNameDataHelper
resource: python/add_other_utilities_to_python.cpp#L80-L92
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/SortListOfFileNameData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`pybind11::list GetSortedListOfFileNameDataHelper( std::vector<FileNameDataCollector::FileNameData>& rFileNameDataList, const std::vector<std::string> & rSortingFlagsOrder )`

# Calls

- [SortListOfFileNameData](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/SortListOfFileNameData.md)