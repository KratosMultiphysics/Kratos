---
type: C++ Function
title: HasSubProperties
resource: includes/properties.h#L477-L480
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/properties/AddSubProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_properties_to_python/HasSubProperties1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/pGetProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_properties/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_materials_input/TestMaterialsInput/_check_results_with_subproperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_properties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool HasSubProperties(const IndexType SubPropertyIndex) const`

# Called by

- [AddSubProperties](../../../functions/includes/properties/AddSubProperties.md)
- [HasSubProperties1](../../../functions/python/add_properties_to_python/HasSubProperties1.md)
- [HasProperties](../../../functions/sources/model_part/ModelPart/HasProperties.md)
- [pGetProperties](../../../functions/sources/model_part/ModelPart/pGetProperties.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/containers/test_properties/KRATOS_TEST_CASE_IN_SUITE.md)
- [_check_results_with_subproperties](../../../functions/tests/test_materials_input/TestMaterialsInput/_check_results_with_subproperties.md)
- [test_model_part_properties](../../../functions/tests/test_model_part/TestModelPart/test_model_part_properties.md)