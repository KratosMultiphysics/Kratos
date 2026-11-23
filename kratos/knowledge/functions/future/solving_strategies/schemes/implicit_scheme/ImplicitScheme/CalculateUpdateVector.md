---
type: C++ Method
title: CalculateUpdateVector
resource: future/solving_strategies/schemes/implicit_scheme.h#L1154-L1170
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Update
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateUpdateVector(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData)`

# Calls

- [pGetVector](../../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)

# Called by

- [Update](../../../../../../functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Update.md)