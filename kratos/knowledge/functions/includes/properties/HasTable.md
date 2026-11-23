---
type: C++ Function
title: HasTable
resource: includes/properties.h#L638-L641
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python/add_properties_to_python/HasTableHelperFunction1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/includes/test_property_accessor/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_properties/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_has_table
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/CheckIfOverwritingTable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool HasTable(const TXVariableType& XVariable, const TYVariableType& YVariable) const`

# Called by

- [HasTableHelperFunction1](../../../functions/python/add_properties_to_python/HasTableHelperFunction1.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/includes/test_property_accessor/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_properties/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_has_table](../../../functions/tests/test_properties/TestProperties/test_has_table.md)
- [CheckIfOverwritingTable](../../../functions/utilities/read_materials_utility/CheckIfOverwritingTable.md)