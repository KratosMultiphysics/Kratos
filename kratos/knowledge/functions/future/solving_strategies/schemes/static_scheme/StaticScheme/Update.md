---
type: C++ Method
title: Update
resource: future/solving_strategies/schemes/static_scheme.h#L219-L248
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/UpdateConstraintsOnlyDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateUpdateVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/GetMoveMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Update(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override`

# Calls

- [pGetVector](../../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)
- [UpdateConstraintsOnlyDofs](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/UpdateConstraintsOnlyDofs.md)
- [CalculateUpdateVector](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateUpdateVector.md)
- [IsFree](../../../../../../functions/includes/dof/Dof/IsFree.md)
- [GetMoveMesh](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/GetMoveMesh.md)