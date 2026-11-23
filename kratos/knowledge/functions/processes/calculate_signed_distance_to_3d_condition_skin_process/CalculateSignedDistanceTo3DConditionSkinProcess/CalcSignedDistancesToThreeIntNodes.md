---
type: C++ Method
title: CalcSignedDistancesToThreeIntNodes
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L950-L963
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToApproxTriangle
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcSignedDistancesToThreeIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element, std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure, std::vector<TetEdgeStruct> IntersectedTetEdges, array_1d<double,4>& ElementalDistances)`

# Calls

- [CalcSignedDistancesToApproxTriangle](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToApproxTriangle.md)