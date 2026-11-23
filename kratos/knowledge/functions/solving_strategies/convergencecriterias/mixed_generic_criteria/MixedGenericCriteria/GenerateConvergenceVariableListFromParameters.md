---
type: C++ Method
title: GenerateConvergenceVariableListFromParameters
resource: solving_strategies/convergencecriterias/mixed_generic_criteria.h#L548-L571
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/MixedGenericCriteria
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static ConvergenceVariableListType GenerateConvergenceVariableListFromParameters(Kratos::Parameters ThisParameters)`

# Calls

- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetDouble](../../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)

# Called by

- [MixedGenericCriteria](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/MixedGenericCriteria.md)