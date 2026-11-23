---
type: C++ Method
title: BreadthFirstSearch
resource: utilities/model_part_graph_utilities.cpp#L180-L203
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeConnectedComponents
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartGraphUtilities::BreadthFirstSearch( const int startVertex, const int color, const DenseVector<ModelPartGraphUtilities::IndexType>& rRowIndices, const DenseVector<ModelPartGraphUtilities::IndexType>& rColIndices, std::unordered_map<ModelPartGraphUtilities::IndexType, int>& visited)`

# Called by

- [ComputeConnectedComponents](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeConnectedComponents.md)