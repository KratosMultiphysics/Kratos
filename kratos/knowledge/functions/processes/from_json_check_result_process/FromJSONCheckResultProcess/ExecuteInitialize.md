---
type: C++ Method
title: ExecuteInitialize
resource: processes/from_json_check_result_process.cpp#L61-L96
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FillVariablesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ComputeRelevantDigits
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FromJSONCheckResultProcess::ExecuteInitialize()`

# Calls

- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [FillVariablesList](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FillVariablesList.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [InitializeDatabases](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [ComputeRelevantDigits](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ComputeRelevantDigits.md)