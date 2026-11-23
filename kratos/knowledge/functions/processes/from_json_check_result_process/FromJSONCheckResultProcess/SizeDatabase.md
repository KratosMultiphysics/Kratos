---
type: C++ Method
title: SizeDatabase
resource: processes/from_json_check_result_process.cpp#L545-L630
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/GetNodeIdentifier
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FromJSONCheckResultProcess::SizeType FromJSONCheckResultProcess::SizeDatabase( const Parameters& rResults, const NodesArrayType& rNodesArray, const ElementsArrayType& rElementsArray )`

# Calls

- [GetNodeIdentifier](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/GetNodeIdentifier.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)

# Called by

- [InitializeDatabases](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases.md)