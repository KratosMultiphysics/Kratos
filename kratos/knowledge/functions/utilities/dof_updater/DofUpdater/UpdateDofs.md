---
type: C++ Method
title: UpdateDofs
resource: utilities/dof_updater.h#L106-L119
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/Update
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void UpdateDofs( DofsArrayType& rDofSet, const SystemVectorType& rDx)`

# Calls

- [IsFree](../../../../functions/includes/dof/Dof/IsFree.md)

# Called by

- [Update](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Update.md)
- [Update](../../../../functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/Update.md)
- [Update](../../../../functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/Update.md)
- [Update](../../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Update.md)
- [Update](../../../../functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/Update.md)