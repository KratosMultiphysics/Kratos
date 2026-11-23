---
type: C++ Method
title: ComputeH
resource: elements/levelset_convection_element_simplex.h#L366-L380
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double ComputeH(BoundedMatrix<double,TNumNodes, TDim>& DN_DX, const double Volume)`

# Called by

- [CalculateLocalSystem](../../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem.md)