---
type: C++ Method
title: GenerateElements
resource: utilities/cutting_utility.h#L659-L919
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateElements (ModelPart& this_model_part, ModelPart& new_model_part, DenseVector<int> Elems_In_Plane, boost::numeric::ublas::compressed_matrix<int>& Coord, const array_1d<double, 3 >& versor, int plane_number)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)