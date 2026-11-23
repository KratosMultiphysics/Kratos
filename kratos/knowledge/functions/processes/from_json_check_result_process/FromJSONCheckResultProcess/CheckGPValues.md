---
type: C++ Method
title: CheckGPValues
resource: processes/from_json_check_result_process.cpp#L456-L540
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/xml_utilities/xml_elements_array/XmlElementsArray/GetElements
    resolved_by: tree-sitter
    confidence: exact
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

`void FromJSONCheckResultProcess::CheckGPValues(IndexType& rCheckCounter)`

# Calls

- [GetElements](../../../../functions/utilities/xml_utilities/xml_elements_array/XmlElementsArray/GetElements.md)
- [GetVariableData](../../../../functions/sources/serializer/Serializer/GetVariableData.md)
- [CheckValues](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckValues.md)
- [FailMessage](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FailMessage.md)

# Called by

- [ExecuteFinalizeSolutionStep](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ExecuteFinalizeSolutionStep.md)