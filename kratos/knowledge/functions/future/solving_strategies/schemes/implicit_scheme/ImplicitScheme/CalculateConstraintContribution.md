---
type: C++ Method
title: CalculateConstraintContribution
resource: future/solving_strategies/schemes/implicit_scheme.h#L1498-L1513
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool CalculateConstraintContribution( ModelPart::MasterSlaveConstraintType& rConstraint, TLSType& rTLS, const ProcessInfo &rProcessInfo)`

# Called by

- [BuildMasterSlaveConstraints](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints.md)