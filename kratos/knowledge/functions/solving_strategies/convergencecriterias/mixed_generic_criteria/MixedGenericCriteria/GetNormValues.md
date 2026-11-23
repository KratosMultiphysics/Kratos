---
type: C++ Method
title: GetNormValues
resource: solving_strategies/convergencecriterias/mixed_generic_criteria.h#L483-L541
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/FindVarLocalKey
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/CalculateConvergenceNorms
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void GetNormValues( const ModelPart& rModelPart, const DofsArrayType& rDofSet, const TSystemVectorType& rDx, std::vector<int>& rDofsCount, std::vector<TDataType>& rSolutionNormsVector, std::vector<TDataType>& rIncreaseNormsVector) const`

# Calls

- [IsFree](../../../../../functions/includes/dof/Dof/IsFree.md)
- [FindVarLocalKey](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/FindVarLocalKey.md)

# Called by

- [CalculateConvergenceNorms](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/CalculateConvergenceNorms.md)