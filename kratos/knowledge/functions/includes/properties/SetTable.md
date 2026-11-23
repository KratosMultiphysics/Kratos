---
type: C++ Function
title: SetTable
resource: includes/properties.h#L430-L433
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python/add_properties_to_python/SetTableHelperFunction1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadTableBlock
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
  - target: functions/tests/test_properties/TestProperties/test_print
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignTablesToProperty
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetTable(const TXVariableType& XVariable, const TYVariableType& YVariable, TableType const& rThisTable)`

# Called by

- [SetTableHelperFunction1](../../../functions/python/add_properties_to_python/SetTableHelperFunction1.md)
- [ReadTableBlock](../../../functions/sources/model_part_io/ModelPartIO/ReadTableBlock.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/includes/test_property_accessor/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_properties/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_has_table](../../../functions/tests/test_properties/TestProperties/test_has_table.md)
- [test_print](../../../functions/tests/test_properties/TestProperties/test_print.md)
- [AssignTablesToProperty](../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignTablesToProperty.md)