---
type: C++ Method
title: CalcNodalDistancesOfTetNodes
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L485-L516
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistanceTo3DSkin
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/DistanceFluidStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcNodalDistancesOfTetNodes( ModelPart::ElementsContainerType::iterator& i_fluidElement, BoundedMatrix<unsigned int,6,2> TetEdgeIndexTable)`

# Calls

- [GetIntersectedLeaves](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves.md)
- [CalcNodalDistanceTo3DSkin](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistanceTo3DSkin.md)

# Called by

- [DistanceFluidStructure](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/DistanceFluidStructure.md)