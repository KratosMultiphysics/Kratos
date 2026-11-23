---
type: C++ Method
title: CalculateLocalSystem
resource: elements/levelset_convection_element_simplex_algebraic_stabilization.cpp#L73-L246
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/outer_prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::CalculateLocalSystem( MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [outer_prod](../../../../functions/includes/amatrix_interface/outer_prod.md)