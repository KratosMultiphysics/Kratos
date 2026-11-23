---
type: C++ Method
title: FirstLoop
resource: utilities/cutting_utility.h#L367-L504
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/cutting_utility/CuttingUtility/GenerateCut
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FirstLoop(ModelPart& this_model_part, boost::numeric::ublas::compressed_matrix<int>& Coord, const array_1d<double, 3 >& versor, const array_1d<double, 3 >& Xp, int number_of_triangles, DenseVector<int>& Elems_In_Plane, double tolerance)//`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [GenerateCut](../../../../functions/utilities/cutting_utility/CuttingUtility/GenerateCut.md)