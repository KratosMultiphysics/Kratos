---
type: C++ Method
title: CheckGeometricalPolynomialDegree
resource: utilities/specifications_utilities.cpp#L1298-L1387
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  called_by:
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

`int SpecificationsUtilities::CheckGeometricalPolynomialDegree(const ModelPart& rModelPart)`

# Calls

- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)

# Called by

- [test_specifications_utilities_elements](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements.md)
- [test_specifications_utilities_conditions](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions.md)
- [test_specifications_utilities_elements_dependencies](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies.md)
- [test_specifications_utilities_conditions_dependencies](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies.md)