---
type: Python Method
title: _ClearObsoleteRestartFiles
resource: python_scripts/restart_utility.py#L196-L207
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/restart_utility/RestartUtility/SaveRestart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _ClearObsoleteRestartFiles(self):`

# Calls

- [_GetRestartFiles](../../../../functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles.md)
- [__GetFolderPathSave](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave.md)
- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)

# Called by

- [SaveRestart](../../../../functions/python_scripts/restart_utility/RestartUtility/SaveRestart.md)