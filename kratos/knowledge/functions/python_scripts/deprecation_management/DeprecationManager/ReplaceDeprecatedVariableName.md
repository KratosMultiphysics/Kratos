---
type: Python Method
title: ReplaceDeprecatedVariableName
resource: python_scripts/deprecation_management.py#L37-L54
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/tikz_output_process/TikZOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/unv_output_process/UnvOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtk_output_process/VtkOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ReplaceDeprecatedVariableName(parameters, old_variable_name, new_variable_name):`

# Calls

- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [IsInt](../../../../functions/sources/kratos_parameters/Parameters/IsInt.md)
- [SetInt](../../../../functions/sources/kratos_parameters/Parameters/SetInt.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [IsDouble](../../../../functions/sources/kratos_parameters/Parameters/IsDouble.md)
- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [IsString](../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [IsBool](../../../../functions/sources/kratos_parameters/Parameters/IsBool.md)
- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [IsStringArray](../../../../functions/sources/kratos_parameters/Parameters/IsStringArray.md)
- [SetStringArray](../../../../functions/sources/kratos_parameters/Parameters/SetStringArray.md)
- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)

# Called by

- [TranslateLegacyVariablesAccordingToCurrentStandard](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)
- [TranslateLegacyVariablesAccordingToCurrentStandard](../../../../functions/python_scripts/tikz_output_process/TikZOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)
- [TranslateLegacyVariablesAccordingToCurrentStandard](../../../../functions/python_scripts/unv_output_process/UnvOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)
- [TranslateLegacyVariablesAccordingToCurrentStandard](../../../../functions/python_scripts/vtk_output_process/VtkOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)