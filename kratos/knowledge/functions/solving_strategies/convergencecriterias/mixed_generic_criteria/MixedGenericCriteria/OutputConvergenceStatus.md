---
type: C++ Method
title: OutputConvergenceStatus
resource: solving_strategies/convergencecriterias/mixed_generic_criteria.h#L352-L375
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/PostCriteria
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void OutputConvergenceStatus( const std::tuple<std::vector<TDataType>, std::vector<TDataType>>& rConvergenceNorms)`

# Called by

- [PostCriteria](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/PostCriteria.md)