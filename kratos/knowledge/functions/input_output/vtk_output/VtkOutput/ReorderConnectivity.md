---
type: C++ Method
title: ReorderConnectivity
resource: input_output/vtk_output.cpp#L1048-L1094
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteConnectivity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`VtkOutput::GeometryType::Pointer VtkOutput::ReorderConnectivity(GeometryType::Pointer& pGeometry) const`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [WriteConnectivity](../../../../functions/input_output/vtk_output/VtkOutput/WriteConnectivity.md)