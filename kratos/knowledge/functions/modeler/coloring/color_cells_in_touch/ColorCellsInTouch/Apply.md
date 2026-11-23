---
type: C++ Method
title: Apply
resource: modeler/coloring/color_cells_in_touch.cpp#L28-L57
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/ApplyColorToCellsInTouchWithGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorCellsInTouch::Apply() const`

# Calls

- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ApplyColorToCellsInTouchWithGeometry](../../../../../functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/ApplyColorToCellsInTouchWithGeometry.md)