---
type: C++ Method
title: CheckValidityOfModelPartsForOperations
resource: utilities/model_part_operation_utilities.cpp#L293-L364
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodesPointerSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodePointersForEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_CheckValidityOfModelPartsForOperations
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ModelPartOperationUtilities::CheckValidityOfModelPartsForOperations( const ModelPart& rMainModelPart, const std::vector<ModelPart const*>& rCheckModelParts, const bool ThrowError)`

# Calls

- [FillNodesPointerSet](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodesPointerSet.md)
- [FillNodePointersForEntities](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodePointersForEntities.md)
- [CheckNodes](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckNodes.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)

# Called by

- [test_CheckValidityOfModelPartsForOperations](../../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_CheckValidityOfModelPartsForOperations.md)