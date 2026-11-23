---
type: C++ Method
title: CalcSignedDistancesToFourIntNodes
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L968-L1008
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/FindIndexNodesOfTriangle2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToApproxTriangle
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistanceTo3DSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcSignedDistancesToFourIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element, std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure, std::vector<TetEdgeStruct> IntersectedTetEdges, array_1d<double,4>& ElementalDistances)`

# Calls

- [FindIndexNodesOfTriangle2](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/FindIndexNodesOfTriangle2.md)
- [CalcSignedDistancesToApproxTriangle](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToApproxTriangle.md)

# Called by

- [CalcNodalDistanceTo3DSkin](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistanceTo3DSkin.md)