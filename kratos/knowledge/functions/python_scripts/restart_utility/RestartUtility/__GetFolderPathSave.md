---
type: Python Method
title: __GetFolderPathSave
resource: python_scripts/restart_utility.py#L257-L261
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/restart_utility/RestartUtility/SaveRestart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetFileNameSave
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __GetFolderPathSave(self):`

# Called by

- [SaveRestart](../../../../functions/python_scripts/restart_utility/RestartUtility/SaveRestart.md)
- [CreateOutputFolder](../../../../functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder.md)
- [_GetRestartFiles](../../../../functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles.md)
- [_ClearObsoleteRestartFiles](../../../../functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles.md)
- [__GetFileNameSave](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetFileNameSave.md)