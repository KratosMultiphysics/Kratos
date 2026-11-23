---
type: C++ Function
title: SplitStringByDelimiter
resource: utilities/string_utilities.cpp#L145-L159
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/utilities/distributed_model_part_initializer/RecursiveCreateModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/distributed_model_part_initializer/CreateSubModelPartHierarchy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry/Registry/RemoveItem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry/Registry/HasItem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateMultipleTypes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> SplitStringByDelimiter( const std::string& rString, const char Delimiter )`

# Called by

- [RecursiveCreateModelParts](../../../functions/mpi/utilities/distributed_model_part_initializer/RecursiveCreateModelParts.md)
- [CreateSubModelPartHierarchy](../../../functions/mpi/utilities/distributed_model_part_initializer/CreateSubModelPartHierarchy.md)
- [RemoveItem](../../../functions/sources/registry/Registry/RemoveItem.md)
- [HasItem](../../../functions/sources/registry/Registry/HasItem.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [GenerateMultipleTypes](../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateMultipleTypes.md)
- [InitializeParser](../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser.md)