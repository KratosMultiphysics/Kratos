---
type: C++ Method
title: ComputeDiscontinuousInterpolation
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L531-L684
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/MappingPressureToStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeDiscontinuousInterpolation( const Node& pNode, Geometry< Node >& geom, const array_1d<double,4>& distances, array_1d<double,4>& Npos, array_1d<double,4>& Nneg)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [MappingPressureToStructure](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/MappingPressureToStructure.md)