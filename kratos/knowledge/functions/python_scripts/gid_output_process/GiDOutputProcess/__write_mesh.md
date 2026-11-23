---
type: Python Method
title: __write_mesh
resource: python_scripts/gid_output_process.py#L491-L503
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/InitializeMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/FinalizeMesh
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteBeforeSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __write_mesh(self, label):`

# Calls

- [InitializeMesh](../../../../functions/includes/gid_io/InitializeMesh.md)
- [FinalizeMesh](../../../../functions/includes/gid_io/FinalizeMesh.md)

# Called by

- [ExecuteBeforeSolutionLoop](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteBeforeSolutionLoop.md)
- [PrintOutput](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/PrintOutput.md)