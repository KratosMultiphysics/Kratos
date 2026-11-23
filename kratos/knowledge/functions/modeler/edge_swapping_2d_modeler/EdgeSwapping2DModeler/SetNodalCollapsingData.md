---
type: C++ Method
title: SetNodalCollapsingData
resource: modeler/edge_swapping_2d_modeler.h#L469-L520
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Distance2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetCollapsingData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetNodalCollapsingData(NodeType& rThisNode, ModelPart& rThisModelPart)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Distance2](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Distance2.md)

# Called by

- [SetCollapsingData](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetCollapsingData.md)