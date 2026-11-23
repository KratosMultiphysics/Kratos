---
type: C++ Method
title: ReadMaterialsUtility
resource: utilities/read_materials_utility.cpp#L56-L84
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ReadMaterialsUtility::ReadMaterialsUtility( Parameters Params, Model& rModel ) : mrModel(rModel)`

# Calls

- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock.md)