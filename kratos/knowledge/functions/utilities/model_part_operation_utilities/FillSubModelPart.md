---
type: C++ Function
title: FillSubModelPart
resource: utilities/model_part_operation_utilities.h#L92-L119
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillOutputSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void FillSubModelPart( ModelPart& rOutputSubModelPart, const ModelPart& rMainModelPart, const std::vector<ModelPart const*>& rOperands, const bool AddNeighbourEntities)`

# Calls

- [FullName](../../../functions/includes/model_part/FullName.md)
- [FillOutputSubModelPart](../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillOutputSubModelPart.md)