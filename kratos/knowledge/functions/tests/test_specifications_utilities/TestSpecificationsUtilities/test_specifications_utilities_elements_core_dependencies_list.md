---
type: Python Method
title: test_specifications_utilities_elements_core_dependencies_list
resource: tests/test_specifications_utilities.py#L185-L208
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/GetNodalSolutionStepDataSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_specifications_utilities_elements_core_dependencies_list(self):`

# Calls

- [GetNodalSolutionStepDataSize](../../../../functions/includes/model_part/GetNodalSolutionStepDataSize.md)
- [AddMissingVariablesFromEntitiesList](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromEntitiesList.md)
- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)
- [AddMissingDofsFromEntitiesList](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromEntitiesList.md)