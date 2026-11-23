---
type: C++ Method
title: ComputeConnectedComponents
resource: utilities/model_part_graph_utilities.cpp#L64-L99
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/BreadthFirstSearch
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::pair<ModelPartGraphUtilities::IndexType, DenseVector<double>> ModelPartGraphUtilities::ComputeConnectedComponents( const ModelPart::NodesContainerType& rNodes, const DenseVector<ModelPartGraphUtilities::IndexType>& rRowIndices, const DenseVector<ModelPartGraphUtilities::IndexType>& rColIndices )`

# Calls

- [BreadthFirstSearch](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/BreadthFirstSearch.md)

# Called by

- [test_graph_and_connected_components_non_consecutive](../../../../functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive.md)