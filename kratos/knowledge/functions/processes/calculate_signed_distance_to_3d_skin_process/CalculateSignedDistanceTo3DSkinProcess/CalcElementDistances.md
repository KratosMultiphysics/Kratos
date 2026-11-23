---
type: C++ Method
title: CalcElementDistances
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L866-L893
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcDistanceTo3DSkin
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/DistanceFluidStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcElementDistances( ModelPart::ElementsContainerType::iterator& i_fluidElement, BoundedMatrix<unsigned int,6,2> TetEdgeIndexTable )`

# Calls

- [GetIntersectedLeaves](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves.md)
- [CalcDistanceTo3DSkin](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcDistanceTo3DSkin.md)

# Called by

- [DistanceFluidStructure](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/DistanceFluidStructure.md)