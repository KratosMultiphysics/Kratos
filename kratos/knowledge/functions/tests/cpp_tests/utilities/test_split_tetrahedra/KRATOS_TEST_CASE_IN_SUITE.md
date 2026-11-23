---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_split_tetrahedra.cpp#L88-L108
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraSplitMode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra/TetrahedraSplit/Split_Tetrahedra
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_split_tetrahedra/SearchExternalFaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_split_tetrahedra/PrintEdgesAndExternalFaces
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(TetrahedraSplitModes, KratosCoreFastSuite)`

# Calls

- [TetrahedraSplitMode](../../../../../functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraSplitMode.md)
- [Split_Tetrahedra](../../../../../functions/utilities/split_tetrahedra/TetrahedraSplit/Split_Tetrahedra.md)
- [SearchExternalFaces](../../../../../functions/tests/cpp_tests/utilities/test_split_tetrahedra/SearchExternalFaces.md)
- [PrintEdgesAndExternalFaces](../../../../../functions/tests/cpp_tests/utilities/test_split_tetrahedra/PrintEdgesAndExternalFaces.md)