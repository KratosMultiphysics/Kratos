---
type: C++ Method
title: AddMissingVariablesFromSpecifications
resource: utilities/specifications_utilities.cpp#L158-L244
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/SetSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SpecificationsUtilities::AddMissingVariablesFromSpecifications( ModelPart& rModelPart, const Parameters SpecificationsParameters, const std::string EntityName )`

# Calls

- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [SetSolutionStepVariablesList](../../../../functions/includes/node/SetSolutionStepVariablesList.md)

# Called by

- [AddMissingVariables](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariables.md)
- [AddMissingVariablesFromEntitiesList](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromEntitiesList.md)