---
type: C++ Function
title: PreComputeExplicitConstraintMassAndInertia
resource: utilities/constraint_utilities.cpp#L296-L406
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PreComputeExplicitConstraintMassAndInertia( ModelPart& rModelPart, const std::string& DofDisplacementVariableName, const std::string& MassVariableName, const std::string& DofRotationVariableName, const std::string& InertiaVariableName )`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)