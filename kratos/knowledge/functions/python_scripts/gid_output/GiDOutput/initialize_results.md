---
type: Python Method
title: initialize_results
resource: python_scripts/gid_output.py#L168-L188
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/gid_output/GiDOutput/get_out_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/_write_mesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/_initialize_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/write_step_to_list
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def initialize_results(self, model_part):`

# Calls

- [get_out_model_part](../../../../functions/python_scripts/gid_output/GiDOutput/get_out_model_part.md)
- [_write_mesh](../../../../functions/python_scripts/gid_output/GiDOutput/_write_mesh.md)
- [_initialize_results](../../../../functions/python_scripts/gid_output/GiDOutput/_initialize_results.md)
- [write_step_to_list](../../../../functions/python_scripts/gid_output/GiDOutput/write_step_to_list.md)