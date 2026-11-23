---
type: Python Method
title: test_graph_and_connected_components_non_consecutive
resource: tests/test_graph_utilities.py#L7-L109
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/Free
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeCSRGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeConnectedComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/exception/Exception/where
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeConnectedComponentsWithActiveNodesCheck
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ApplyMinimalScalarFixity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_graph_and_connected_components_non_consecutive(self):`

# Calls

- [Free](../../../../functions/includes/node/Free.md)
- [ComputeCSRGraph](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeCSRGraph.md)
- [ComputeConnectedComponents](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeConnectedComponents.md)
- [where](../../../../functions/sources/exception/Exception/where.md)
- [ComputeConnectedComponentsWithActiveNodesCheck](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeConnectedComponentsWithActiveNodesCheck.md)
- [Fix](../../../../functions/includes/node/Fix.md)
- [ApplyMinimalScalarFixity](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ApplyMinimalScalarFixity.md)