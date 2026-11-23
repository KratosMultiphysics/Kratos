---
type: C++ Function
title: SetCommonColumn
resource: utilities/result_dabatase.h#L545-L548
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataTXT
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataJSON
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_result_database/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetCommonColumn(const Vector& rCommonColumn)`

# Called by

- [ReadDataTXT](../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataTXT.md)
- [ReadDataJSON](../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataJSON.md)
- [InitializeDatabases](../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_result_database/KRATOS_TEST_CASE_IN_SUITE.md)