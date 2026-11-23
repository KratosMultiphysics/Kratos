---
type: C++ Method
title: GetBrepSurfaceParameters
resource: input_output/cad_json_output.cpp#L43-L134
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetParameters
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CadJsonOutput::GetBrepSurfaceParameters( const typename ModelPart::GeometryIterator& rGeometryIterator, Parameters& rBrepsParameters, IndexType EchoLevel)`

# Calls

- [AddInt](../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [AddBool](../../../../functions/sources/kratos_parameters/Parameters/AddBool.md)
- [AddEmptyArray](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyArray.md)
- [KnotsU](../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [KnotsV](../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [fill](../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [AddString](../../../../functions/sources/kratos_parameters/Parameters/AddString.md)
- [GetBoundaryLoopParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters.md)

# Called by

- [GetParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetParameters.md)