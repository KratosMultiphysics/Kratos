---
type: C++ Method
title: CalculateLocalSystem
resource: elements/levelset_convection_element_simplex.h#L122-L255
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/ComputeH
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/outer_prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ComputeH](../../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/ComputeH.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [outer_prod](../../../../functions/includes/amatrix_interface/outer_prod.md)