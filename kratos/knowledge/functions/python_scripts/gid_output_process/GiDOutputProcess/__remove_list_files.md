---
type: Python Method
title: __remove_list_files
resource: python_scripts/gid_output_process.py#L661-L674
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/exists
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__restart_list_files
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __remove_list_files(self):`

# Calls

- [exists](../../../../functions/sources/kratos_filesystem/exists.md)
- [remove](../../../../functions/sources/kratos_filesystem/remove.md)

# Called by

- [__restart_list_files](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__restart_list_files.md)