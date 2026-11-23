---
type: C++ Method
title: GetVariableData
resource: sources/serializer.cpp#L31-L34
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/kratos_components/pGet
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataTXT
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataJSON
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckNodeValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckNodeHistoricalValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckGPValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FillDatabase
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_result_database/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`const VariableData* Serializer::GetVariableData(std::string const & VariableName)`

# Calls

- [pGet](../../../../functions/includes/kratos_components/pGet.md)

# Called by

- [ExecuteInitializeSolutionStep](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ExecuteInitializeSolutionStep.md)
- [ReadDataTXT](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataTXT.md)
- [ReadDataJSON](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataJSON.md)
- [CheckNodeValues](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckNodeValues.md)
- [CheckNodeHistoricalValues](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckNodeHistoricalValues.md)
- [CheckGPValues](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckGPValues.md)
- [FillDatabase](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FillDatabase.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_result_database/KRATOS_TEST_CASE_IN_SUITE.md)