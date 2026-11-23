---
type: C++ Function
title: ReplaceAllSubstrings
resource: utilities/string_utilities.cpp#L164-L177
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateTemplatedTypes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/BasicGenericFunctionUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/FunctionBody
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string ReplaceAllSubstrings( const std::string& rInputString, const std::string& rStringToBeReplaced, const std::string& rStringToReplace )`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [GenerateTemplatedTypes](../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateTemplatedTypes.md)
- [BasicGenericFunctionUtility](../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/BasicGenericFunctionUtility.md)
- [FunctionBody](../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/FunctionBody.md)
- [InitializeParser](../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser.md)