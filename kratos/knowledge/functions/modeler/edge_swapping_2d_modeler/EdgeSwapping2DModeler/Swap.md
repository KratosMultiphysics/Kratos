---
type: C++ Method
title: Swap
resource: modeler/edge_swapping_2d_modeler.h#L532-L552
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Remesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Swap(Element& rElement1, Element& rElement2, int Edge1, int Edge2)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [Remesh](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Remesh.md)