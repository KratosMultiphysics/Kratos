---
type: C++ Method
title: GetNodeIdentifier
resource: processes/from_json_check_result_process.cpp#L811-L814
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/SizeDatabase
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FillDatabase
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string FromJSONCheckResultProcess::GetNodeIdentifier(NodeType& rNode)`

# Called by

- [SizeDatabase](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/SizeDatabase.md)
- [FillDatabase](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FillDatabase.md)