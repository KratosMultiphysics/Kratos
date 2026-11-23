---
type: C++ Method
title: UpdateDerivatives
resource: solving_strategies/schemes/residual_based_bdf_scheme.h#L377-L397
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/Update
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void UpdateDerivatives( ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb )`

# Called by

- [Update](../../../../../functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/Update.md)