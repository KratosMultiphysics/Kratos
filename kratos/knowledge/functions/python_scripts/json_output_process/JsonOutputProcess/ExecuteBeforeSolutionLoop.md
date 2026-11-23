---
type: Python Method
title: ExecuteBeforeSolutionLoop
resource: python_scripts/json_output_process.py#L95-L229
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/__check_flag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/__get_node_identifier
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/json_utilities/write_external_json
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteBeforeSolutionLoop(self):`

# Calls

- [__check_flag](../../../../functions/python_scripts/json_output_process/JsonOutputProcess/__check_flag.md)
- [__get_node_identifier](../../../../functions/python_scripts/json_output_process/JsonOutputProcess/__get_node_identifier.md)
- [write_external_json](../../../../functions/python_scripts/json_utilities/write_external_json.md)