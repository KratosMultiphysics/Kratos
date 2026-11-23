---
type: C++ Method
title: InitializeSolutionStep
resource: solving_strategies/schemes/residual_based_bdf_scheme.h#L236-L258
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/time_discretization/BDF/ComputeAndSaveBDFCoefficients
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeSolutionStep( ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [ComputeAndSaveBDFCoefficients](../../../../../functions/utilities/time_discretization/BDF/ComputeAndSaveBDFCoefficients.md)