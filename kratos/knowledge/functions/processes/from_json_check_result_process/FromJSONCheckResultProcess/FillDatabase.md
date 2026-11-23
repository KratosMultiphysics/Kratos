---
type: C++ Method
title: FillDatabase
resource: processes/from_json_check_result_process.cpp#L635-L806
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/serializer/Serializer/GetVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/GetNodeIdentifier
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FromJSONCheckResultProcess::FillDatabase( const Parameters& rResults, const NodesArrayType& rNodesArray, const ElementsArrayType& rElementsArray, const SizeType NumberOfGP )`

# Calls

- [GetVariableData](../../../../functions/sources/serializer/Serializer/GetVariableData.md)
- [GetNodeIdentifier](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/GetNodeIdentifier.md)