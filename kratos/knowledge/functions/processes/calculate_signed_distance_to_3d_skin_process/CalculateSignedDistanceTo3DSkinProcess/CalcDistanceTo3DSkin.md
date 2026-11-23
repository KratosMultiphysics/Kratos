---
type: C++ Method
title: CalcDistanceTo3DSkin
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1145-L1189
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/FillIntNodesContainer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AvoidZeroDistances
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcElementDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcDistanceTo3DSkin( std::vector<TetEdgeStruct>& IntersectedTetEdges, ModelPart::ElementsContainerType::iterator& i_fluid_Element, unsigned int NumberIntersectionsOnTetCorner )`

# Calls

- [FillIntNodesContainer](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/FillIntNodesContainer.md)
- [CalcSignedDistancesToMoreThanThreeIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes.md)
- [AvoidZeroDistances](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AvoidZeroDistances.md)

# Called by

- [CalcElementDistances](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcElementDistances.md)