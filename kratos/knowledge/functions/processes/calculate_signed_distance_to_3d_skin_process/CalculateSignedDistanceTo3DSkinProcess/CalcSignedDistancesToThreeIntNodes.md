---
type: C++ Method
title: CalcSignedDistancesToThreeIntNodes
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1275-L1309
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PointDistanceToPlane
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcSignedDistancesToThreeIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_Element, std::vector<IntersectionNodeStruct>& NodesOfApproximatedStructure, array_1d<double,4>& ElementalDistances )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointDistanceToPlane](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PointDistanceToPlane.md)