---
type: C++ Method
title: FindVarLocalKey
resource: solving_strategies/convergencecriterias/mixed_generic_criteria.h#L415-L429
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variable_data/IsComponent
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/GetNormValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool FindVarLocalKey( typename DofsArrayType::const_iterator itDof, int& rVarLocalKey) const`

# Calls

- [GetVariable](../../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [IsComponent](../../../../../functions/containers/variable_data/IsComponent.md)

# Called by

- [GetNormValues](../../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/GetNormValues.md)