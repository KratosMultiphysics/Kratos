---
type: C++ Function
title: PreComputeExplicitConstraintConstribution
resource: utilities/constraint_utilities.cpp#L182-L291
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

`void PreComputeExplicitConstraintConstribution( ModelPart& rModelPart, const std::vector<std::string>& rDofVariableNames, const std::vector<std::string>& rResidualDofVariableNames )`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)