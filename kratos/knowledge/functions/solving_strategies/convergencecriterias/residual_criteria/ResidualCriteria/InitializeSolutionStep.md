---
type: C++ Method
title: InitializeSolutionStep
resource: solving_strategies/convergencecriterias/residual_criteria.h#L234-L251
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/CalculateResidualNorm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeSolutionStep( ModelPart& rModelPart, DofsArrayType& rDofSet, const TSystemMatrixType& rA, const TSystemVectorType& rDx, const TSystemVectorType& rb ) override`

# Calls

- [CalculateResidualNorm](../../../../../functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/CalculateResidualNorm.md)