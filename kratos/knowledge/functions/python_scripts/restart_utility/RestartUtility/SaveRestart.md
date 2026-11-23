---
type: Python Method
title: SaveRestart
resource: python_scripts/restart_utility.py#L119-L155
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetPrettyTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetFileNameSave
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/save_restart_process/SaveRestartProcess/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/__execute_restart_utility_save
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def SaveRestart(self):`

# Calls

- [__GetPrettyTime](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetPrettyTime.md)
- [CreateOutputFolder](../../../../functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder.md)
- [__GetFolderPathSave](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave.md)
- [__GetFileNameSave](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetFileNameSave.md)
- [_ClearObsoleteRestartFiles](../../../../functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles.md)

# Called by

- [PrintOutput](../../../../functions/python_scripts/save_restart_process/SaveRestartProcess/PrintOutput.md)
- [__execute_restart_utility_save](../../../../functions/tests/test_restart/TestRestart/__execute_restart_utility_save.md)