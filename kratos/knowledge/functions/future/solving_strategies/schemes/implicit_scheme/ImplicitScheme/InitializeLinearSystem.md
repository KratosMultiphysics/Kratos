---
type: C++ Method
title: InitializeLinearSystem
resource: future/solving_strategies/schemes/implicit_scheme.h#L1362-L1391
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpDofArrays
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpSystemIds
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeLinearSystem(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData)`

# Calls

- [SetUpDofArrays](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpDofArrays.md)
- [SetUpSystemIds](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpSystemIds.md)

# Called by

- [Predict](../../../../../../functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Predict.md)