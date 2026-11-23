---
type: C++ Method
title: GetIntersectedLeaves
resource: spatial_containers/octree_binary.h#L1001-L1088
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/HasChildren
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMinPointNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMaxPointNormalized
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistancesOfTetNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcElementDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/FindIntersectedSkinObjects
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/IdentifyNearEntitiesAndCheckEntityForIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GetIntersectedLeaves( typename cell_type::pointer_type pObject, std::vector<cell_type*>& rLeaves, const double ToleranceCoefficient = 0.001 // 0.1% of the min size )`

# Calls

- [CalcKeyNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized.md)
- [pop_back](../../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)
- [HasChildren](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/HasChildren.md)
- [GetMinPointNormalized](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMinPointNormalized.md)
- [GetMaxPointNormalized](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMaxPointNormalized.md)

# Called by

- [CalcNodalDistancesOfTetNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistancesOfTetNodes.md)
- [CalcElementDistances](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcElementDistances.md)
- [FindIntersectedSkinObjects](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/FindIntersectedSkinObjects.md)
- [IdentifyNearEntitiesAndCheckEntityForIntersection](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/IdentifyNearEntitiesAndCheckEntityForIntersection.md)