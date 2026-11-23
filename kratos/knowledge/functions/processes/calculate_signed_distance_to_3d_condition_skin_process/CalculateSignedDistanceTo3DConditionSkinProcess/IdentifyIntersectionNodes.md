---
type: C++ Method
title: IdentifyIntersectionNodes
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L521-L631
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IdentifyIntersectionNodes( ModelPart::ElementsContainerType::iterator& i_fluidElement, unsigned int i_tetEdge, std::vector<OctreeType::cell_type*>& leaves, std::vector<TetEdgeStruct>& IntersectedTetEdges, unsigned int& NumberIntersectionsOnTetCorner, BoundedMatrix<unsigned int,6,2> TetEdgeIndexTable, int& intersection_counter)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)