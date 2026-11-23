---
type: C++ Method
title: SensitivityBuilder
resource: utilities/sensitivity_builder.cpp#L319-L329
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sensitivity_builder/ReplaceDeprecatedNameIfExists
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/pGetSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/GetVariableLists
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SensitivityBuilder::SensitivityBuilder( Parameters Settings, ModelPart& rModelPart, AdjointResponseFunction::Pointer pResponseFunction) : SensitivityBuilder( Settings, rModelPart, pResponseFunction, Kratos::make_unique<SensitivityBuilderScheme>())`

# Calls

- [ReplaceDeprecatedNameIfExists](../../../../functions/utilities/sensitivity_builder/ReplaceDeprecatedNameIfExists.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [pGetSubModelPart](../../../../functions/sources/model_part/ModelPart/pGetSubModelPart.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetVariableLists](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/GetVariableLists.md)
- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)