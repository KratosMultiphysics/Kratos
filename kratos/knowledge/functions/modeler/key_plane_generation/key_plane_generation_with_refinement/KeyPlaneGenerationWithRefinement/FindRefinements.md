---
type: C++ Method
title: FindRefinements
resource: modeler/key_plane_generation/key_plane_generation_with_refinement.cpp#L96-L118
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void KeyPlaneGenerationWithRefinement::FindRefinements(std::vector<BoundingBoxType>& rRefinementAreas, std::vector<array_1d<double,3>>& rVoxelSizes,const array_1d<double,3>& rGlobalVoxelSize)`

# Calls

- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [Generate](../../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/Generate.md)