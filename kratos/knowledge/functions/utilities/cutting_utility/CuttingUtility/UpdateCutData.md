---
type: C++ Method
title: UpdateCutData
resource: utilities/cutting_utility.h#L930-L949
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/GetNodalSolutionStepDataSize
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/write_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nodal_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nonhistorical_nodal_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nodal_flags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__write_elemental_conditional_flags
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void UpdateCutData( ModelPart& new_model_part, ModelPart& old_model_part)`

# Calls

- [GetNodalSolutionStepDataSize](../../../../functions/includes/model_part/GetNodalSolutionStepDataSize.md)

# Called by

- [write_results](../../../../functions/python_scripts/gid_output/GiDOutput/write_results.md)
- [__write_nodal_results](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nodal_results.md)
- [__write_nonhistorical_nodal_results](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nonhistorical_nodal_results.md)
- [__write_nodal_flags](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nodal_flags.md)
- [__write_elemental_conditional_flags](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_elemental_conditional_flags.md)