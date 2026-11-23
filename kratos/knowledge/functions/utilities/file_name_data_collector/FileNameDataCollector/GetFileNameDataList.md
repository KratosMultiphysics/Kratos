---
type: C++ Method
title: GetFileNameDataList
resource: utilities/file_name_data_collector.cpp#L164-L183
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/ListDirectory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/filename
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/RetrieveFileNameData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetSortedFileNamesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<FileNameDataCollector::FileNameData> FileNameDataCollector::GetFileNameDataList() const`

# Calls

- [ListDirectory](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/ListDirectory.md)
- [filename](../../../../functions/sources/kratos_filesystem/filename.md)
- [RetrieveFileNameData](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/RetrieveFileNameData.md)

# Called by

- [_GetRestartFiles](../../../../functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles.md)
- [GetSortedFileNamesList](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetSortedFileNamesList.md)