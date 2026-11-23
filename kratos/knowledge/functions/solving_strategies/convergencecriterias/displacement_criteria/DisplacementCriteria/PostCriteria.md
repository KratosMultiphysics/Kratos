---
type: C++ Method
title: PostCriteria
resource: solving_strategies/convergencecriterias/displacement_criteria.h#L179-L225
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/CalculateFinalCorrectionNorm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/CalculateReferenceNorm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool PostCriteria( ModelPart& rModelPart, DofsArrayType& rDofSet, const TSystemMatrixType& rA, const TSystemVectorType& rDx, const TSystemVectorType& rb ) override`

# Calls

- [CalculateFinalCorrectionNorm](../../../../../functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/CalculateFinalCorrectionNorm.md)
- [CalculateReferenceNorm](../../../../../functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/CalculateReferenceNorm.md)