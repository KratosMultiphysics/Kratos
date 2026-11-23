---
type: C++ Method
title: GetDofsListFromSpecifications
resource: utilities/specifications_utilities.cpp#L377-L390
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromElementsSpecifications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromConditionsSpecifications
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_GetDofsListFromSpecifications
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> SpecificationsUtilities::GetDofsListFromSpecifications(const ModelPart& rModelPart)`

# Calls

- [GetDofsListFromElementsSpecifications](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromElementsSpecifications.md)
- [GetDofsListFromConditionsSpecifications](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromConditionsSpecifications.md)

# Called by

- [test_specifications_utilities_GetDofsListFromSpecifications](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_GetDofsListFromSpecifications.md)