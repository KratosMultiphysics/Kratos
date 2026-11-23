---
type: C++ Method
title: SetLoadState
resource: sources/serializer.cpp#L38-L41
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/serializer/Serializer/SeekBegin
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_serializer_to_python/SerializerLoadFromBeginning
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_serializer/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_unique_ptr_aliasing/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Serializer::SetLoadState()`

# Calls

- [SeekBegin](../../../../functions/sources/serializer/Serializer/SeekBegin.md)

# Called by

- [SerializerLoadFromBeginning](../../../../functions/python/add_serializer_to_python/SerializerLoadFromBeginning.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_serializer/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_unique_ptr_aliasing/KRATOS_TEST_CASE_IN_SUITE.md)