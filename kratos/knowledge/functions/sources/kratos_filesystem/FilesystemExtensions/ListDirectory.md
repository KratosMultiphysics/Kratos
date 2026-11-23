---
type: C++ Method
title: ListDirectory
resource: sources/kratos_filesystem.cpp#L115-L122
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::filesystem::path> FilesystemExtensions::ListDirectory(const std::filesystem::path& rPath)`

# Called by

- [GetFileNameDataList](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList.md)