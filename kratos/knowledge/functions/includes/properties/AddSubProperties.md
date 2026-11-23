---
type: C++ Function
title: AddSubProperties
resource: includes/properties.h#L466-L470
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/properties/HasSubProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/containers/test_properties/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_properties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_sub_properties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddSubProperties(Properties::Pointer pNewSubProperty)`

# Calls

- [HasSubProperties](../../../functions/includes/properties/HasSubProperties.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/containers/test_properties/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_model_part_properties](../../../functions/tests/test_model_part/TestModelPart/test_model_part_properties.md)
- [test_model_part_sub_properties](../../../functions/tests/test_model_part/TestModelPart/test_model_part_sub_properties.md)