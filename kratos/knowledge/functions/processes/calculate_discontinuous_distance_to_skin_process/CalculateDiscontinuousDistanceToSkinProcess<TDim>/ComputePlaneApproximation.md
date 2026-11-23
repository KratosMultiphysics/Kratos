---
type: C++ Method
title: ComputePlaneApproximation
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L551-L565
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputePlaneApproximation( const Element& rElement1, const std::vector< array_1d<double,3> >& rPointsCoord, array_1d<double,3>& rPlaneBasePointCoords, array_1d<double,3>& rPlaneNormal)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)