---
type: C++ Method
title: RecursiveAddEntities
resource: utilities/model_part_combination_utilities.cpp#L246-L271
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/FastTransferBetweenModelPartsProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartCombinationUtilities::RecursiveAddEntities( ModelPart& rDestinationModelPart, ModelPart& rOriginModelPart )`

# Calls

- [FastTransferBetweenModelPartsProcess](../../../../functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/FastTransferBetweenModelPartsProcess.md)
- [NumberOfSubModelParts](../../../../functions/includes/model_part/NumberOfSubModelParts.md)