---
type: C++ Method
title: IsEquivalentTo
resource: sources/kratos_parameters.cpp#L1240-L1287
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_recursive_validation_4_levels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Parameters::IsEquivalentTo(Parameters& rParameters)`

# Calls

- [key](../../../../functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_recursive_validation_4_levels](../../../../functions/tests/test_kratos_parameters/TestParameters/test_recursive_validation_4_levels.md)
- [test_specifications_utilities_elements](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements.md)
- [test_specifications_utilities_conditions](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions.md)
- [test_specifications_utilities_elements_dependencies](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies.md)
- [test_specifications_utilities_conditions_dependencies](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies.md)