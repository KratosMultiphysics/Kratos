---
type: C++ Function
title: for_pure_c11
resource: utilities/parallel_utilities.h#L486-L512
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/check_scalar_base_process/aux_object_cpp_callback/f
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/exception/Exception/what
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_parallel_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void for_pure_c11(TUnaryFunction &&f)`

# Calls

- [f](../../../functions/python_scripts/check_scalar_base_process/aux_object_cpp_callback/f.md)
- [what](../../../functions/sources/exception/Exception/what.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_parallel_utilities/KRATOS_TEST_CASE_IN_SUITE.md)