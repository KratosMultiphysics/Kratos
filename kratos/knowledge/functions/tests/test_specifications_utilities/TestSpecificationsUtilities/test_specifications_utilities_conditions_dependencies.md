---
type: Python Method
title: test_specifications_utilities_conditions_dependencies
resource: tests/test_specifications_utilities.py#L242-L276
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/GetNodalSolutionStepDataSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineTimeIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineFramework
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineSymmetricLHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DeterminePositiveDefiniteLHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineIfCompatibleGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineIfRequiresTimeIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/CheckCompatibleConstitutiveLaws
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/CheckGeometricalPolynomialDegree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/GetDocumention
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsEquivalentTo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_specifications_utilities_conditions_dependencies(self):`

# Calls

- [GetNodalSolutionStepDataSize](../../../../functions/includes/model_part/GetNodalSolutionStepDataSize.md)
- [AddMissingVariables](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariables.md)
- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)
- [AddMissingDofs](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofs.md)
- [DetermineTimeIntegration](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineTimeIntegration.md)
- [DetermineFramework](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineFramework.md)
- [DetermineSymmetricLHS](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineSymmetricLHS.md)
- [DeterminePositiveDefiniteLHS](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DeterminePositiveDefiniteLHS.md)
- [DetermineIfCompatibleGeometries](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineIfCompatibleGeometries.md)
- [DetermineIfRequiresTimeIntegration](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineIfRequiresTimeIntegration.md)
- [CheckCompatibleConstitutiveLaws](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/CheckCompatibleConstitutiveLaws.md)
- [CheckGeometricalPolynomialDegree](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/CheckGeometricalPolynomialDegree.md)
- [GetDocumention](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/GetDocumention.md)
- [IsEquivalentTo](../../../../functions/sources/kratos_parameters/Parameters/IsEquivalentTo.md)