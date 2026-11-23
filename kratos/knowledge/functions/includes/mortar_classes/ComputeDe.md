---
type: C++ Function
title: ComputeDe
resource: includes/mortar_classes.h#L1538-L1553
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAeComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAeComponents
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GeometryMatrixType ComputeDe( const TArray& N1, const double detJ ) const`

# Called by

- [CalculateAeComponents](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAeComponents.md)
- [CalculateDeltaAeComponents](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAeComponents.md)