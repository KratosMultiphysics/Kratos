---
type: C++ Method
title: CheckOriginModelPartsAndAssignRoot
resource: modeler/combine_model_part_modeler.cpp#L83-L106
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::CheckOriginModelPartsAndAssignRoot()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/SetupModelPart.md)