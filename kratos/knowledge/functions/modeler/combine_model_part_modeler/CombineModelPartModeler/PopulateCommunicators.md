---
type: C++ Method
title: PopulateCommunicators
resource: modeler/combine_model_part_modeler.cpp#L199-L222
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateLocalMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::PopulateCommunicators()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [PopulateLocalMesh](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateLocalMesh.md)
- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/SetupModelPart.md)