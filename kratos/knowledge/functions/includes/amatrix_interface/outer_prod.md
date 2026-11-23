---
type: C++ Function
title: outer_prod
resource: includes/amatrix_interface.h#L582-L588
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex_algebraic_stabilization/LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAeComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAeComponents
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`AMatrix::VectorOuterProductExpression<TExpression1Type, TExpression2Type> outer_prod( AMatrix::MatrixExpression<TExpression1Type, TCategory1> const& First, AMatrix::MatrixExpression<TExpression2Type, TCategory2> const& Second)`

# Called by

- [CalculateLocalSystem](../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem.md)
- [CalculateLocalSystem](../../../functions/elements/levelset_convection_element_simplex_algebraic_stabilization/LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>/CalculateLocalSystem.md)
- [CalculateAeComponents](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAeComponents.md)
- [CalculateDeltaAeComponents](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAeComponents.md)