---
type: C++ Method
title: ApplyMinimalScalarFixity
resource: utilities/model_part_graph_utilities.cpp#L147-L178
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<ModelPartGraphUtilities::IndexType> ModelPartGraphUtilities::ApplyMinimalScalarFixity( ModelPart::NodesContainerType& rNodes, const Variable<double>& rVar, const DenseVector<double>& colors, const ModelPartGraphUtilities::IndexType ncolors )`

# Calls

- [Fix](../../../../functions/includes/node/Fix.md)

# Called by

- [test_graph_and_connected_components_non_consecutive](../../../../functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive.md)