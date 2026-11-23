---
type: C++ Method
title: PostCriteria
resource: solving_strategies/convergencecriterias/mixed_generic_criteria.h#L166-L187
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/CalculateConvergenceNorms
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/OutputConvergenceStatus
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/CheckConvergence
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool PostCriteria( ModelPart& rModelPart, DofsArrayType& rDofSet, const TSystemMatrixType& A, const TSystemVectorType& Dx, const TSystemVectorType& b) override`

# Calls

- [CalculateConvergenceNorms](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/CalculateConvergenceNorms.md)
- [OutputConvergenceStatus](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/OutputConvergenceStatus.md)
- [CheckConvergence](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/CheckConvergence.md)