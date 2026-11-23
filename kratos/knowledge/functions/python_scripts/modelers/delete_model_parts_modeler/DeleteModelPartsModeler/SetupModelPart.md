---
type: Python Method
title: SetupModelPart
resource: python_scripts/modelers/delete_model_parts_modeler.py#L21-L54
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/__GetSubModelPartsNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/DeleteModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def SetupModelPart(self):`

# Calls

- [GetStringArray](../../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [HasModelPart](../../../../../functions/containers/model/Model/HasModelPart.md)
- [__GetSubModelPartsNames](../../../../../functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/__GetSubModelPartsNames.md)
- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [DeleteModelPart](../../../../../functions/containers/model/Model/DeleteModelPart.md)