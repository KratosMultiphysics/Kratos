---
type: C++ Method
title: QuadraticErrorFunctionPoint
resource: utilities/qef_utility.cpp#L28-L65
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/UpdateATAATB
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/ComputeQuadraticErrorFunctionPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/FirstEnd
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/SecondEnd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`array_1d<double, 3> QuadraticErrorFunction::QuadraticErrorFunctionPoint ( const GeometryType& rVoxel, const GeometryArrayType& rTriangles )`

# Calls

- [UpdateATAATB](../../../../functions/utilities/qef_utility/QuadraticErrorFunction/UpdateATAATB.md)
- [ComputeQuadraticErrorFunctionPoint](../../../../functions/utilities/qef_utility/QuadraticErrorFunction/ComputeQuadraticErrorFunctionPoint.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [FirstEnd](../../../../functions/utilities/qef_utility/QuadraticErrorFunction/FirstEnd.md)
- [SecondEnd](../../../../functions/utilities/qef_utility/QuadraticErrorFunction/SecondEnd.md)