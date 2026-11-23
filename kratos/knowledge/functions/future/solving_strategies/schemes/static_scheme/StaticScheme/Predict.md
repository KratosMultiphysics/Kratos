---
type: C++ Method
title: Predict
resource: future/solving_strategies/schemes/static_scheme.h#L175-L217
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeLinearSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/builder/Builder/CalculateSolutionVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/GetMoveMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Predict(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override`

# Calls

- [InitializeLinearSystem](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeLinearSystem.md)
- [CalculateSolutionVector](../../../../../../functions/future/solving_strategies/builders/builder/Builder/CalculateSolutionVector.md)
- [GetMoveMesh](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/GetMoveMesh.md)