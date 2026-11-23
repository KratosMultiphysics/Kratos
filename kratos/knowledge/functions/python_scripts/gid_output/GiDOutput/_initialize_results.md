---
type: Python Method
title: _initialize_results
resource: python_scripts/gid_output.py#L156-L157
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/InitializeResults
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

`def _initialize_results(self, label, model_part):`

# Calls

- [InitializeResults](../../../../functions/includes/gid_io/InitializeResults.md)

# Called by

- [initialize_results](../../../../functions/python_scripts/gid_output/GiDOutput/initialize_results.md)
- [write_results](../../../../functions/python_scripts/gid_output/GiDOutput/write_results.md)