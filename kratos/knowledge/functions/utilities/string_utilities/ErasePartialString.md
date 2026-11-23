---
type: C++ Function
title: ErasePartialString
resource: utilities/string_utilities.cpp#L84-L101
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataTXT
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string ErasePartialString( const std::string& rMainString, const std::string& rToErase )`

# Called by

- [IdentifyDataTXT](../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataTXT.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE.md)