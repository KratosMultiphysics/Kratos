---
type: C++ Method
title: Predict
resource: solving_strategies/schemes/residual_based_bdf_custom_scheme.h#L286-L328
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ComputePredictComponent
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Predict( ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b ) override`

# Calls

- [ComputePredictComponent](../../../../../functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ComputePredictComponent.md)