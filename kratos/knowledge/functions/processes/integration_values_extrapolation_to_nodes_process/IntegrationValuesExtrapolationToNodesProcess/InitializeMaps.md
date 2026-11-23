---
type: C++ Method
title: InitializeMaps
resource: processes/integration_values_extrapolation_to_nodes_process.cpp#L295-L368
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/ExecuteBeforeSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationValuesExtrapolationToNodesProcess::InitializeMaps()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [ExecuteBeforeSolutionLoop](../../../../functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/ExecuteBeforeSolutionLoop.md)