---
type: Python Method
title: __write_elemental_conditional_flags
resource: python_scripts/gid_output_process.py#L551-L559
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/PrintFlagsOnGaussPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/UpdateCutData
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

`def __write_elemental_conditional_flags(self, label):`

# Calls

- [PrintFlagsOnGaussPoints](../../../../functions/includes/gid_io/PrintFlagsOnGaussPoints.md)
- [UpdateCutData](../../../../functions/utilities/cutting_utility/CuttingUtility/UpdateCutData.md)

# Called by

- [ExecuteBeforeSolutionLoop](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteBeforeSolutionLoop.md)
- [PrintOutput](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/PrintOutput.md)