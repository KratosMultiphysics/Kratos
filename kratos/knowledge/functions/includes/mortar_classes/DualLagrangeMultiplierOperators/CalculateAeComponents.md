---
type: C++ Method
title: CalculateAeComponents
resource: includes/mortar_classes.h#L1488-L1499
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/mortar_classes/ComputeDe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/outer_prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAeComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CalculateAe
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateAeComponents( KinematicVariables& rKinematicVariables, const double rIntegrationWeight )`

# Calls

- [ComputeDe](../../../../functions/includes/mortar_classes/ComputeDe.md)
- [outer_prod](../../../../functions/includes/amatrix_interface/outer_prod.md)

# Called by

- [CalculateDeltaAeComponents](../../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAeComponents.md)
- [CalculateAe](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CalculateAe.md)