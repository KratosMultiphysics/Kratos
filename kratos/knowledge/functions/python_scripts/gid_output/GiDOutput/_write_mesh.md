---
type: Python Method
title: _write_mesh
resource: python_scripts/gid_output.py#L151-L154
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
  - target: functions/python_scripts/gid_output/GiDOutput/initialize_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/write_results
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _write_mesh(self, label, model_part):`

# Calls

- [InitializeMesh](../../../../functions/includes/gid_io/InitializeMesh.md)
- [FinalizeMesh](../../../../functions/includes/gid_io/FinalizeMesh.md)

# Called by

- [initialize_results](../../../../functions/python_scripts/gid_output/GiDOutput/initialize_results.md)
- [write_results](../../../../functions/python_scripts/gid_output/GiDOutput/write_results.md)