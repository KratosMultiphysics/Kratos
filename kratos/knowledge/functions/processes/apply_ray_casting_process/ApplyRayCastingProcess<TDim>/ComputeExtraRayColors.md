---
type: C++ Method
title: ComputeExtraRayColors
resource: processes/apply_ray_casting_process.cpp#L213-L268
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/CorrectExtraRayOrigin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/GetRayIntersections
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/DistancePositionInSpace
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyRayCastingProcess<TDim>::ComputeExtraRayColors( const array_1d<double,3> &rCoords, array_1d<double,TDim> &rDistances)`

# Calls

- [CorrectExtraRayOrigin](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/CorrectExtraRayOrigin.md)
- [GetRayIntersections](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/GetRayIntersections.md)

# Called by

- [DistancePositionInSpace](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/DistancePositionInSpace.md)