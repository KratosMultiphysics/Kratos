---
type: C++ Method
title: CreateCommunicators
resource: modeler/combine_model_part_modeler.cpp#L166-L194
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/SetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::CreateCommunicators()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetNumberOfColors](../../../../functions/sources/communicator/Communicator/SetNumberOfColors.md)
- [GetNumberOfColors](../../../../functions/sources/communicator/Communicator/GetNumberOfColors.md)
- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)