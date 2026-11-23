---
type: C++ Method
title: CalculateConvergenceNorms
resource: solving_strategies/convergencecriterias/mixed_generic_criteria.h#L308-L345
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/GetNormValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/PostCriteria
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::tuple<std::vector<TDataType>, std::vector<TDataType>> CalculateConvergenceNorms( const ModelPart& rModelPart, const DofsArrayType& rDofSet, const TSystemVectorType& rDx)`

# Calls

- [GetNormValues](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/GetNormValues.md)

# Called by

- [PostCriteria](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/PostCriteria.md)