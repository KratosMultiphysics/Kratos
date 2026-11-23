---
type: Python Method
title: CreateOutputFolder
resource: python_scripts/restart_utility.py#L166-L168
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/restart_utility/RestartUtility/SaveRestart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/save_restart_process/SaveRestartProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CreateOutputFolder(self):`

# Calls

- [MPISafeCreateDirectories](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)
- [__GetFolderPathSave](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave.md)

# Called by

- [SaveRestart](../../../../functions/python_scripts/restart_utility/RestartUtility/SaveRestart.md)
- [__init__](../../../../functions/python_scripts/save_restart_process/SaveRestartProcess/__init__.md)