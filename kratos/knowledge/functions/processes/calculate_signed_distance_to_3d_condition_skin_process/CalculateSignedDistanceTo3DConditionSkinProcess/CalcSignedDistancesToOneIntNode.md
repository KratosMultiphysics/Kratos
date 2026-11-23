---
type: C++ Method
title: CalcSignedDistancesToOneIntNode
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L854-L890
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalcSignedDistancesToOneIntNode( ModelPart::ElementsContainerType::iterator& i_fluid_element, std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure, array_1d<double,4>& ElementalDistances)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)