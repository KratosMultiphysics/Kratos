---
type: C++ Method
title: AddMissingDofsFromSpecifications
resource: utilities/specifications_utilities.cpp#L341-L372
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SpecificationsUtilities::AddMissingDofsFromSpecifications( ModelPart& rModelPart, const Parameters SpecificationsParameters, const std::string EntityName )`

# Calls

- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)

# Called by

- [AddMissingDofs](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofs.md)
- [AddMissingDofsFromEntitiesList](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromEntitiesList.md)