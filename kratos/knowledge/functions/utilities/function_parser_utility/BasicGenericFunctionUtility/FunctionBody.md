---
type: C++ Method
title: FunctionBody
resource: utilities/function_parser_utility.cpp#L94-L102
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/ReplaceAllSubstrings
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_function_parser_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_GenericFunctionUtility1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_GenericFunctionUtility2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtility
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string BasicGenericFunctionUtility::FunctionBody()`

# Calls

- [ReplaceAllSubstrings](../../../../functions/utilities/string_utilities/ReplaceAllSubstrings.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_function_parser_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_GenericFunctionUtility1](../../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_GenericFunctionUtility1.md)
- [test_GenericFunctionUtility2](../../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_GenericFunctionUtility2.md)
- [test_ApplyFunctionToNodesUtility](../../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtility.md)