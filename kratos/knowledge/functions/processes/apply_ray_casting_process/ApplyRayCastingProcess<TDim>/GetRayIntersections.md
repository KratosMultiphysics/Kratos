---
type: C++ Method
title: GetRayIntersections
resource: processes/apply_ray_casting_process.cpp#L335-L385
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/DistancePositionInSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/ComputeExtraRayColors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyRayCastingProcess<TDim>::GetRayIntersections( const double* ray, const unsigned int direction, std::vector<std::pair<double,Element::GeometryType*> >& rIntersections)`

# Calls

- [CalcKeyNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized.md)

# Called by

- [DistancePositionInSpace](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/DistancePositionInSpace.md)
- [ComputeExtraRayColors](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/ComputeExtraRayColors.md)