---
type: Python Method
title: __GetFileNameSave
resource: python_scripts/restart_utility.py#L268-L271
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/restart_utility/RestartUtility/SaveRestart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __GetFileNameSave(self, file_label):`

# Calls

- [__GetFolderPathSave](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave.md)

# Called by

- [SaveRestart](../../../../functions/python_scripts/restart_utility/RestartUtility/SaveRestart.md)