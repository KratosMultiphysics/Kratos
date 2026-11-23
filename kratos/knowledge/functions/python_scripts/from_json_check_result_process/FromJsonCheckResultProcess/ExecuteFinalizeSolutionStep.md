---
type: Python Method
title: ExecuteFinalizeSolutionStep
resource: python_scripts/from_json_check_result_process.py#L42-L51
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/from_json_check_result_process/IsCorrectResult
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/GetErrorMessage
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteFinalizeSolutionStep(self):`

# Calls

- [IsCorrectResult](../../../../functions/processes/from_json_check_result_process/IsCorrectResult.md)
- [GetErrorMessage](../../../../functions/processes/from_json_check_result_process/GetErrorMessage.md)