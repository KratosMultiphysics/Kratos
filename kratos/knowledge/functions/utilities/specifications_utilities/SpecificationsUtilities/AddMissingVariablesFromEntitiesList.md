---
type: C++ Method
title: AddMissingVariablesFromEntitiesList
resource: utilities/specifications_utilities.cpp#L126-L153
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromSpecifications
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_core_dependencies_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies_list
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SpecificationsUtilities::AddMissingVariablesFromEntitiesList( ModelPart& rModelPart, const Parameters EntitiesList )`

# Calls

- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [AddMissingVariablesFromSpecifications](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromSpecifications.md)

# Called by

- [test_specifications_utilities_elements_list](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_list.md)
- [test_specifications_utilities_conditions_list](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_list.md)
- [test_specifications_utilities_elements_core_dependencies_list](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_core_dependencies_list.md)
- [test_specifications_utilities_elements_dependencies_list](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies_list.md)
- [test_specifications_utilities_conditions_dependencies_list](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies_list.md)