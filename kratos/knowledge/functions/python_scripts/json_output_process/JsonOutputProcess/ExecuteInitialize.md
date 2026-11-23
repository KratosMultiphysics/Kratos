---
type: Python Method
title: ExecuteInitialize
resource: python_scripts/json_output_process.py#L62-L93
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/TotalProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/__generate_variable_list_from_input
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteInitialize(self):`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [TotalProcesses](../../../../functions/sources/communicator/Communicator/TotalProcesses.md)
- [__generate_variable_list_from_input](../../../../functions/python_scripts/json_output_process/JsonOutputProcess/__generate_variable_list_from_input.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)