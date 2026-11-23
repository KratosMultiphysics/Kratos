---
type: Python Method
title: _GetRestartFiles
resource: python_scripts/restart_utility.py#L172-L191
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameData/GetStep
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _GetRestartFiles(self):`

# Calls

- [__GetFolderPathSave](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave.md)
- [FileNameDataCollector](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector.md)
- [GetFileNameDataList](../../../../functions/utilities/file_name_data_collector/FileNameDataCollector/GetFileNameDataList.md)
- [GetStep](../../../../functions/utilities/file_name_data_collector/FileNameData/GetStep.md)

# Called by

- [_ClearObsoleteRestartFiles](../../../../functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles.md)