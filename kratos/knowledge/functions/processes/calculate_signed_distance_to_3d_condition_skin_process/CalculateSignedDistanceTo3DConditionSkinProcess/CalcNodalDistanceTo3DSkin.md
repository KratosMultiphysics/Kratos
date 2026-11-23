---
type: C++ Method
title: CalcNodalDistanceTo3DSkin
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L770-L811
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/ComputeApproximationNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToFourIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/AssignDistancesToElements
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistancesOfTetNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcNodalDistanceTo3DSkin(std::vector<TetEdgeStruct>& IntersectedTetEdges, ModelPart::ElementsContainerType::iterator& i_fluid_element, unsigned int NumberIntersectionsOnTetCorner)`

# Calls

- [ComputeApproximationNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/ComputeApproximationNodes.md)
- [CalcSignedDistancesToFourIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToFourIntNodes.md)
- [AssignDistancesToElements](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/AssignDistancesToElements.md)

# Called by

- [CalcNodalDistancesOfTetNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistancesOfTetNodes.md)