---
type: C++ Method
title: ComputeGraph
resource: utilities/model_part_graph_utilities.cpp#L25-L52
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeCSRGraph
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Kratos::unique_ptr<SparseContiguousRowGraph<>> ModelPartGraphUtilities::ComputeGraph(const ModelPart& rModelPart)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [ComputeCSRGraph](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeCSRGraph.md)