---
type: C++ Method
title: GetGeometryConnectivity
resource: input_output/ensight_output.cpp#L2846-L2898
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EnSightOutput::GetGeometryConnectivity( const GeometricalObject& rGeometricalObject, std::vector<std::size_t>& rConnectivity )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [WriteGeometryFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile.md)