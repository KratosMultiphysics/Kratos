---
type: C++ Method
title: RecursiveAddOfModelPartsToList
resource: utilities/model_part_combination_utilities.cpp#L219-L241
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/model_part_combination_utilities/ModelPartCombinationUtilities/CheckSubModelParts
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartCombinationUtilities::RecursiveAddOfModelPartsToList( ModelPart& rModelPart, std::unordered_map<std::string, std::size_t>& rListModelParts )`

# Calls

- [NumberOfSubModelParts](../../../../functions/includes/model_part/NumberOfSubModelParts.md)

# Called by

- [CheckSubModelParts](../../../../functions/utilities/model_part_combination_utilities/ModelPartCombinationUtilities/CheckSubModelParts.md)