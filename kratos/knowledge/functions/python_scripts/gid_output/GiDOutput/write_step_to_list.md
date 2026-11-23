---
type: Python Method
title: write_step_to_list
resource: python_scripts/gid_output.py#L48-L64
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/initialize_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/write_results
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def write_step_to_list(self, step_label=None):`

# Called by

- [initialize_results](../../../../functions/python_scripts/gid_output/GiDOutput/initialize_results.md)
- [write_results](../../../../functions/python_scripts/gid_output/GiDOutput/write_results.md)