---
type: C++ Method
title: DuplicateMesh
resource: modeler/combine_model_part_modeler.cpp#L131-L161
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::DuplicateMesh() const`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/SetupModelPart.md)