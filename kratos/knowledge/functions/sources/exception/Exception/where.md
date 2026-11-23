---
type: C++ Method
title: where
resource: sources/exception.cpp#L85-L90
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/exception/Exception/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_exceptions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`const CodeLocation Exception::where() const`

# Called by

- [PrintData](../../../../functions/sources/exception/Exception/PrintData.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_exceptions/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_graph_and_connected_components_non_consecutive](../../../../functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive.md)