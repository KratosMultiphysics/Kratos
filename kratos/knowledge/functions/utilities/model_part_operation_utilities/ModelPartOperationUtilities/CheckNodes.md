---
type: C++ Method
title: CheckNodes
resource: utilities/model_part_operation_utilities.cpp#L239-L257
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckValidityOfModelPartsForOperations
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::CheckNodes( std::vector<IndexType>& rNodeIdsWithIssues, const ModelPart::NodesContainerType& rCheckNodes, const std::set<ModelPart::NodeType const*>& rMainNodes)`

# Called by

- [CheckValidityOfModelPartsForOperations](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckValidityOfModelPartsForOperations.md)