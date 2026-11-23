---
type: C++ Method
title: CopyPropertiesModeler
resource: modeler/copy_properties_modeler.cpp#L25-L32
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CopyPropertiesModeler::CopyPropertiesModeler( Model& rModel, Parameters ModelerParameters) : Modeler(rModel, ModelerParameters) , mpModel(&rModel)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)
- [AddString](../../../../functions/sources/kratos_parameters/Parameters/AddString.md)