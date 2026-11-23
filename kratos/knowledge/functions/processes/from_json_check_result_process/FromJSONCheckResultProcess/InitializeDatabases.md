---
type: C++ Method
title: InitializeDatabases
resource: processes/from_json_check_result_process.cpp#L174-L215
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_elements_array/XmlElementsArray/GetElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/SizeDatabase
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/result_dabatase/SetCommonColumn
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FromJSONCheckResultProcess::InitializeDatabases()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetElements](../../../../functions/utilities/xml_utilities/xml_elements_array/XmlElementsArray/GetElements.md)
- [SizeDatabase](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/SizeDatabase.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [SetCommonColumn](../../../../functions/utilities/result_dabatase/SetCommonColumn.md)

# Called by

- [ExecuteInitialize](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ExecuteInitialize.md)