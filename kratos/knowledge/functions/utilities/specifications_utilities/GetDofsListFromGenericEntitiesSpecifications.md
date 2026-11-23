---
type: C++ Function
title: GetDofsListFromGenericEntitiesSpecifications
resource: utilities/specifications_utilities.cpp#L31-L60
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromElementsSpecifications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromConditionsSpecifications
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> GetDofsListFromGenericEntitiesSpecifications(const TContainerType& rContainer)`

# Calls

- [GetStringArray](../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)

# Called by

- [GetDofsListFromElementsSpecifications](../../../functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromElementsSpecifications.md)
- [GetDofsListFromConditionsSpecifications](../../../functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromConditionsSpecifications.md)