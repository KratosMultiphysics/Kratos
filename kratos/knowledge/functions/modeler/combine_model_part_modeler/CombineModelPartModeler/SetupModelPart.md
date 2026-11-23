---
type: C++ Method
title: SetupModelPart
resource: modeler/combine_model_part_modeler.cpp#L40-L63
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
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CheckOriginModelPartsAndAssignRoot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateCommunicators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CreateSubModelParts
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::SetupModelPart()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)
- [CheckOriginModelPartsAndAssignRoot](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CheckOriginModelPartsAndAssignRoot.md)
- [DuplicateMesh](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateMesh.md)
- [PopulateCommunicators](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateCommunicators.md)
- [CreateSubModelParts](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CreateSubModelParts.md)