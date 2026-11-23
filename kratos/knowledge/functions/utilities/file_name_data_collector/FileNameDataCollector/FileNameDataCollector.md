---
type: C++ Method
title: FileNameDataCollector
resource: utilities/file_name_data_collector.cpp#L27-L116
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetPatternFlagStrings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/FindAndReplace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/parent_path
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/filename
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/IsFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/GetPatternFlagString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_GetFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_file_name_data_collector/check_file_names
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_SortListOfFileNameData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_GetSortedFileNamesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FileNameDataCollector::FileNameDataCollector( const ModelPart& rModelPart, const std::string& rFileNamePattern, const std::unordered_map<std::string, std::string>& rFlagFormatMap) : mrModelPart(rModelPart)`

# Calls

- [GetPatternFlagStrings](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetPatternFlagStrings.md)
- [FindAndReplace](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/FindAndReplace.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [parent_path](../../../../functions/sources/kratos_filesystem/parent_path.md)
- [filename](../../../../functions/sources/kratos_filesystem/filename.md)
- [emplace_back](../../../../functions/containers/pointer_vector/emplace_back.md)
- [IsFlag](../../../../functions/utilities/file_name_data_collector/IsFlag.md)
- [GetPatternFlagString](../../../../functions/utilities/file_name_data_collector/GetPatternFlagString.md)

# Called by

- [_GetRestartFiles](../../../../functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles.md)
- [test_GetFileName](../../../../functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_GetFileName.md)
- [check_file_names](../../../../functions/tests/test_file_name_data_collector/check_file_names.md)
- [test_SortListOfFileNameData](../../../../functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_SortListOfFileNameData.md)
- [test_GetSortedFileNamesList](../../../../functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_GetSortedFileNamesList.md)