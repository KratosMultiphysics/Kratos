---
type: C++ Method
title: CalcSignedDistancesToTwoIntNodes
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L895-L945
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcSignedDistancesToTwoIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element, std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure, array_1d<double,4>& ElementalDistances)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointDistanceToLineSegment3D](../../../../functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D.md)