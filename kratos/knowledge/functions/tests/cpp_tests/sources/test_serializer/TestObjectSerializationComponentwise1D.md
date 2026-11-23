---
type: C++ Function
title: TestObjectSerializationComponentwise1D
resource: tests/cpp_tests/sources/test_serializer.cpp#L52-L60
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/sources/test_serializer/SaveAndLoadObjects
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_serializer/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TestObjectSerializationComponentwise1D(const TObjectType& rObjectToBeSaved, TObjectType& rObjectToBeLoaded)`

# Calls

- [SaveAndLoadObjects](../../../../../functions/tests/cpp_tests/sources/test_serializer/SaveAndLoadObjects.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/sources/test_serializer/KRATOS_TEST_CASE_IN_SUITE.md)