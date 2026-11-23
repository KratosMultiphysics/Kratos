---
type: C++ Method
title: FillOutputSubModelPart
resource: utilities/model_part_operation_utilities.cpp#L180-L237
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/model_part_operation_utilities/FillSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::FillOutputSubModelPart( ModelPart& rOutputSubModelPart, ModelPart& rMainModelPart, std::vector<ModelPart::NodeType*>& rOutputNodes, std::vector<ModelPart::ConditionType*>& rOutputConditions, std::vector<ModelPart::ElementType*>& rOutputElements)`

# Calls

- [FullName](../../../../functions/includes/model_part/FullName.md)

# Called by

- [FillSubModelPart](../../../../functions/utilities/model_part_operation_utilities/FillSubModelPart.md)