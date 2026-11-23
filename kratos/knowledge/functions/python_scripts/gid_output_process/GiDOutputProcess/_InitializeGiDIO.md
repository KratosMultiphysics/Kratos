---
type: Python Method
title: _InitializeGiDIO
resource: python_scripts/gid_output_process.py#L339-L362
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__get_gidpost_flag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/GidIO
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _InitializeGiDIO(self,gidpost_flags,param):`

# Calls

- [__get_gidpost_flag](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__get_gidpost_flag.md)
- [GidIO](../../../../functions/includes/gid_io/GidIO.md)

# Called by

- [ExecuteInitialize](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize.md)