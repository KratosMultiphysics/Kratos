---
type: C++ Method
title: CalculateDeltaAeComponents
resource: includes/mortar_classes.h#L1703-L1725
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAeComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/ComputeDe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/outer_prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDeltaAeComponents( KinematicVariablesType& rKinematicVariables, DerivativeDataType& rDerivativeData, const double rIntegrationWeight )`

# Calls

- [CalculateAeComponents](../../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAeComponents.md)
- [ComputeDe](../../../../functions/includes/mortar_classes/ComputeDe.md)
- [outer_prod](../../../../functions/includes/amatrix_interface/outer_prod.md)