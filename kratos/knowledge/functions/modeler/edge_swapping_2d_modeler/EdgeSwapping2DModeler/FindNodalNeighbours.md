---
type: C++ Method
title: FindNodalNeighbours
resource: modeler/edge_swapping_2d_modeler.h#L600-L627
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/CollapseNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetSwappingData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindNodalNeighbours(ModelPart& rThisModelPart)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CollapseNodes](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/CollapseNodes.md)
- [SetSwappingData](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetSwappingData.md)