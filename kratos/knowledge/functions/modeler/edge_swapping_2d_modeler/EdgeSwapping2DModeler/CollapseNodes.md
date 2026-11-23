---
type: C++ Method
title: CollapseNodes
resource: modeler/edge_swapping_2d_modeler.h#L375-L443
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNodalNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetCollapsingData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CollapseNodes(ModelPart& rThisModelPart)`

# Calls

- [FindNodalNeighbours](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNodalNeighbours.md)
- [SetCollapsingData](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetCollapsingData.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)