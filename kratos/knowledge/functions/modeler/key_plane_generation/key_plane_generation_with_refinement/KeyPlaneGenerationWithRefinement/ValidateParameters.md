---
type: C++ Method
title: ValidateParameters
resource: modeler/key_plane_generation/key_plane_generation_with_refinement.cpp#L22-L38
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/GetDefaultRefinementZoneParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void KeyPlaneGenerationWithRefinement::ValidateParameters()`

# Calls

- [ValidateAndAssignDefaults](../../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetDefaultRefinementZoneParameters](../../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/GetDefaultRefinementZoneParameters.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)