---
type: Python Method
title: __finalize_results
resource: python_scripts/gid_output_process.py#L561-L566
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/FinalizeResults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteBeforeSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __finalize_results(self):`

# Calls

- [FinalizeResults](../../../../functions/includes/gid_io/FinalizeResults.md)

# Called by

- [ExecuteBeforeSolutionLoop](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteBeforeSolutionLoop.md)
- [PrintOutput](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/PrintOutput.md)
- [ExecuteFinalize](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteFinalize.md)