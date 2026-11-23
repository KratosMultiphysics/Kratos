---
type: C++ Function
title: GetSensitivityVariableName
resource: utilities/sensitivity_utilities.h#L85-L110
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/variable_data/IsComponent
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesGetSensitivityVariableName/test_GetSensitivityVariableName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/SensitivityVariables/SensitivityVariables
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::string GetSensitivityVariableName( const TVariableType& rVariable)`

# Calls

- [IsComponent](../../../functions/containers/variable_data/IsComponent.md)

# Called by

- [test_GetSensitivityVariableName](../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesGetSensitivityVariableName/test_GetSensitivityVariableName.md)
- [SensitivityVariables](../../../functions/utilities/sensitivity_builder/SensitivityVariables/SensitivityVariables.md)