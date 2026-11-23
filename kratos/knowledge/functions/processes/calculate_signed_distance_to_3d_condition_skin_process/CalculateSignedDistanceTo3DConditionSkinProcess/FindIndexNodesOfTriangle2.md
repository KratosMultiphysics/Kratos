---
type: C++ Method
title: FindIndexNodesOfTriangle2
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1013-L1051
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToFourIntNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindIndexNodesOfTriangle2(std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure, array_1d<unsigned int,3>& IndexNodes_T2)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CalcSignedDistancesToFourIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToFourIntNodes.md)