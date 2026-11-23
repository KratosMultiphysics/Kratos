---
type: C++ Method
title: GetDofsListFromConditionsSpecifications
resource: utilities/specifications_utilities.cpp#L404-L408
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/specifications_utilities/GetDofsListFromGenericEntitiesSpecifications
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromSpecifications
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> SpecificationsUtilities::GetDofsListFromConditionsSpecifications(const ModelPart& rModelPart)`

# Calls

- [GetDofsListFromGenericEntitiesSpecifications](../../../../functions/utilities/specifications_utilities/GetDofsListFromGenericEntitiesSpecifications.md)

# Called by

- [GetDofsListFromSpecifications](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/GetDofsListFromSpecifications.md)