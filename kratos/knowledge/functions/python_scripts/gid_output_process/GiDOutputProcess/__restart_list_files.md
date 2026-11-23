---
type: Python Method
title: __restart_list_files
resource: python_scripts/gid_output_process.py#L597-L659
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_list_files
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_SetCurrentTimeParameters
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __restart_list_files(self,additional_frequencies):`

# Calls

- [__remove_list_files](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_list_files.md)

# Called by

- [_SetCurrentTimeParameters](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_SetCurrentTimeParameters.md)