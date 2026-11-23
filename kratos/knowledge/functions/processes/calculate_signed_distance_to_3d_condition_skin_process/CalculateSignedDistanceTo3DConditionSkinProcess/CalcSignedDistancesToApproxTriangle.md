---
type: C++ Method
title: CalcSignedDistancesToApproxTriangle
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1056-L1132
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToFourIntNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcSignedDistancesToApproxTriangle( ModelPart::ElementsContainerType::iterator& i_fluid_element, std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure, std::vector<TetEdgeStruct> IntersectedTetEdges, array_1d<double,4>& ElementalDistances, array_1d<unsigned int,3> IndexNodes)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CalcSignedDistancesToThreeIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToThreeIntNodes.md)
- [CalcSignedDistancesToFourIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToFourIntNodes.md)