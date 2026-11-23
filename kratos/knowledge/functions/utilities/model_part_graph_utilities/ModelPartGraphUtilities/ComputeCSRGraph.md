---
type: C++ Method
title: ComputeCSRGraph
resource: utilities/model_part_graph_utilities.cpp#L54-L62
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeGraph
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::pair<DenseVector<ModelPartGraphUtilities::IndexType>, DenseVector<ModelPartGraphUtilities::IndexType>> ModelPartGraphUtilities::ComputeCSRGraph(const ModelPart& rModelPart)`

# Calls

- [ComputeGraph](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeGraph.md)

# Called by

- [test_graph_and_connected_components_non_consecutive](../../../../functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive.md)