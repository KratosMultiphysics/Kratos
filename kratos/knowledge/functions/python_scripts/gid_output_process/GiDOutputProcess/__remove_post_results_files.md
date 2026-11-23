---
type: Python Method
title: __remove_post_results_files
resource: python_scripts/gid_output_process.py#L676-L717
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
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_SetCurrentTimeParameters
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __remove_post_results_files(self, step_label):`

# Calls

- [exists](../../../../functions/sources/kratos_filesystem/exists.md)
- [remove](../../../../functions/sources/kratos_filesystem/remove.md)

# Called by

- [_SetCurrentTimeParameters](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_SetCurrentTimeParameters.md)