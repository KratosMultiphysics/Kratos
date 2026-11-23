---
type: C++ Function
title: HasAccessor
resource: includes/properties.h#L412-L415
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python/add_properties_to_python/AddInterfaceToAccessorFold
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/includes/test_property_accessor/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_accessor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_and_set_accessors_utilities/ReadAndSetAccessorsUtilities/ReadAndSetAccessors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool HasAccessor(const TVariableType& rVariable) const`

# Called by

- [AddInterfaceToAccessorFold](../../../functions/python/add_properties_to_python/AddInterfaceToAccessorFold.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/includes/test_property_accessor/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_accessor](../../../functions/tests/test_properties/TestProperties/test_accessor.md)
- [ReadAndSetAccessors](../../../functions/utilities/read_and_set_accessors_utilities/ReadAndSetAccessorsUtilities/ReadAndSetAccessors.md)