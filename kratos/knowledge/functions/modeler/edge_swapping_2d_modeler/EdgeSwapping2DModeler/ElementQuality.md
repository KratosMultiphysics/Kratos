---
type: C++ Method
title: ElementQuality
resource: modeler/edge_swapping_2d_modeler.h#L589-L598
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/SideLenghts2D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/MarkBadQualityElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double ElementQuality(Element& rThisElement)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SideLenghts2D](../../../../functions/utilities/geometry_utilities/SideLenghts2D.md)

# Called by

- [MarkBadQualityElements](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/MarkBadQualityElements.md)