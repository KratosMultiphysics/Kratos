---
type: C++ Method
title: CalculateNodalRayColors
resource: modeler/internals/cartesian_mesh_colors.h#L247-L307
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CalculateColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateNodalRayColors(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, int InsideColor, int OutsideColor)`

# Calls

- [CollapseIntersectionPoints](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints.md)
- [CalculateColor](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CalculateColor.md)

# Called by

- [Execute](../../../../../functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute.md)