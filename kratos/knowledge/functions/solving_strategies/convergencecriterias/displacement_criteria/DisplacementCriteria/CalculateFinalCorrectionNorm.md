---
type: C++ Method
title: CalculateFinalCorrectionNorm
resource: solving_strategies/convergencecriterias/displacement_criteria.h#L462-L497
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/PostCriteria
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`TDataType CalculateFinalCorrectionNorm( SizeType& rDofNum, DofsArrayType& rDofSet, const TSystemVectorType& rDx, ModelPart& rModelPart )`

# Called by

- [PostCriteria](../../../../../functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/PostCriteria.md)