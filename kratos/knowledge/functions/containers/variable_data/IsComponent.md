---
type: C++ Function
title: IsComponent
resource: containers/variable_data.h#L220-L223
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/variable/Variable/Print
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variable/Variable/Info
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list/Add
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list/Has
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/MixedGenericCriteria
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/FindVarLocalKey
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_variables/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_utilities/GetSensitivityVariableName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`constexpr bool IsComponent() const noexcept`

# Called by

- [Print](../../../functions/containers/variable/Variable/Print.md)
- [Info](../../../functions/containers/variable/Variable/Info.md)
- [Add](../../../functions/containers/variables_list/Add.md)
- [Has](../../../functions/containers/variables_list/Has.md)
- [MixedGenericCriteria](../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/MixedGenericCriteria.md)
- [FindVarLocalKey](../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/FindVarLocalKey.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/containers/test_variables/KRATOS_TEST_CASE_IN_SUITE.md)
- [GetSensitivityVariableName](../../../functions/utilities/sensitivity_utilities/GetSensitivityVariableName.md)