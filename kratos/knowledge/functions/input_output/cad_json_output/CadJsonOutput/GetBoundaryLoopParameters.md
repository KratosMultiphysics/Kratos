---
type: C++ Method
title: GetBoundaryLoopParameters
resource: input_output/cad_json_output.cpp#L136-L191
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/HasSameCurveDirection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_curve/BrepCurve/Knots
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CadJsonOutput::GetBoundaryLoopParameters( const BrepCurveOnSurfaceArrayType& rCurveOnSurfaceArray, Parameters& rBoundaryLoop, IndexType EchoLevel)`

# Calls

- [AddEmptyArray](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyArray.md)
- [AddInt](../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [AddBool](../../../../functions/sources/kratos_parameters/Parameters/AddBool.md)
- [HasSameCurveDirection](../../../../functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/HasSameCurveDirection.md)
- [Knots](../../../../functions/geometries/brep_curve/BrepCurve/Knots.md)
- [GetT0](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0.md)
- [GetT1](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1.md)
- [fill](../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)