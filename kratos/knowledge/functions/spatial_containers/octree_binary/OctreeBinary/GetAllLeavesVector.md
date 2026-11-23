---
type: C++ Method
title: GetAllLeavesVector
resource: spatial_containers/octree_binary.h#L469-L485
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/HasChildren
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/PrintGiDMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/PrintGiDResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PrintGiDMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PrintGiDResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/Constrain2To1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/Constrain2To1New
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/PrintGiDMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/PrintGiDMeshNew
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int GetAllLeavesVector(std::vector<cell_type*>& all_leaves) const`

# Calls

- [pop_back](../../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)
- [HasChildren](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/HasChildren.md)

# Called by

- [GenerateNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateNodes.md)
- [CalculateDistance2](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance2.md)
- [CalculateDistance](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance.md)
- [PrintGiDMesh](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/PrintGiDMesh.md)
- [PrintGiDResults](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/PrintGiDResults.md)
- [GenerateNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateNodes.md)
- [CalculateDistance2](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance2.md)
- [CalculateDistance](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance.md)
- [PrintGiDMesh](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PrintGiDMesh.md)
- [PrintGiDResults](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PrintGiDResults.md)
- [Constrain2To1](../../../../functions/spatial_containers/octree_binary/OctreeBinary/Constrain2To1.md)
- [Constrain2To1New](../../../../functions/spatial_containers/octree_binary/OctreeBinary/Constrain2To1New.md)
- [PrintGiDMesh](../../../../functions/spatial_containers/octree_binary/OctreeBinary/PrintGiDMesh.md)
- [PrintGiDMeshNew](../../../../functions/spatial_containers/octree_binary/OctreeBinary/PrintGiDMeshNew.md)