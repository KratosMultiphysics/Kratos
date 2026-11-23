---
type: C++ Method
title: BreadthFirstSearchWithActiveNodesCheck
resource: utilities/model_part_graph_utilities.cpp#L205-L230
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeConnectedComponentsWithActiveNodesCheck
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartGraphUtilities::BreadthFirstSearchWithActiveNodesCheck( const int startVertex, const int color, const DenseVector<ModelPartGraphUtilities::IndexType>& rRowIndices, const DenseVector<ModelPartGraphUtilities::IndexType>& rColIndices, std::unordered_map<ModelPartGraphUtilities::IndexType, int>& visited, const std::unordered_map<ModelPartGraphUtilities::IndexType, bool>& active_nodes)`

# Called by

- [ComputeConnectedComponentsWithActiveNodesCheck](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeConnectedComponentsWithActiveNodesCheck.md)