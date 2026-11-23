---
type: C++ Function
title: ContainsPartialString
resource: utilities/string_utilities.cpp#L106-L123
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataTXT
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ContainsPartialString( const std::string& rMainString, const std::string& rToCheck )`

# Called by

- [AssignScalarInputToEntitiesProcess](../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess.md)
- [IdentifyDataTXT](../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataTXT.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [InitializeParser](../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser.md)