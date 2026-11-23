---
type: C++ Method
title: CheckNodeHistoricalValues
resource: processes/from_json_check_result_process.cpp#L380-L451
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/serializer/Serializer/GetVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FailMessage
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ExecuteFinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FromJSONCheckResultProcess::CheckNodeHistoricalValues(IndexType& rCheckCounter)`

# Calls

- [GetVariableData](../../../../functions/sources/serializer/Serializer/GetVariableData.md)
- [CheckValues](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckValues.md)
- [FailMessage](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FailMessage.md)

# Called by

- [ExecuteFinalizeSolutionStep](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ExecuteFinalizeSolutionStep.md)