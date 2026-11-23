---
type: C++ Method
title: CalcSignedDistancesToMoreThanThreeIntNodes
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1314-L1426
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/EigenVectors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PointDistanceToPlane
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcDistanceTo3DSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcSignedDistancesToMoreThanThreeIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_Element, std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure, array_1d<double,4>& ElementalDistances, std::vector<TetEdgeStruct>& IntersectedTetEdges )`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [EigenVectors](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/EigenVectors.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [PointDistanceToPlane](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PointDistanceToPlane.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CalcDistanceTo3DSkin](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcDistanceTo3DSkin.md)