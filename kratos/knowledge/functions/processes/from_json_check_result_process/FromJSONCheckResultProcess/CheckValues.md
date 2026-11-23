---
type: C++ Method
title: CheckValues
resource: processes/from_json_check_result_process.cpp#L259-L266
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckNodeValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckNodeHistoricalValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckGPValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool FromJSONCheckResultProcess::CheckValues( const double ValueEntity, const double ValueJSON )`

# Called by

- [CheckNodeValues](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckNodeValues.md)
- [CheckNodeHistoricalValues](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckNodeHistoricalValues.md)
- [CheckGPValues](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckGPValues.md)