---
type: Python Function
title: write_external_json
resource: python_scripts/json_utilities.py#L8-L10
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/ExecuteBeforeSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/ExecuteFinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def write_external_json(file_name, data):`

# Called by

- [ExecuteBeforeSolutionLoop](../../../functions/python_scripts/json_output_process/JsonOutputProcess/ExecuteBeforeSolutionLoop.md)
- [ExecuteFinalizeSolutionStep](../../../functions/python_scripts/json_output_process/JsonOutputProcess/ExecuteFinalizeSolutionStep.md)