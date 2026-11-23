---
type: Python Method
title: __write_step_to_list
resource: python_scripts/gid_output_process.py#L568-L595
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteBeforeSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __write_step_to_list(self,step_label=None):`

# Called by

- [ExecuteBeforeSolutionLoop](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteBeforeSolutionLoop.md)
- [PrintOutput](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/PrintOutput.md)