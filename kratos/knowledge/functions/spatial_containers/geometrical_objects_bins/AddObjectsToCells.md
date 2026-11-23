---
type: C++ Function
title: AddObjectsToCells
resource: spatial_containers/geometrical_objects_bins.h#L422-L442
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/GetCellBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddObjectsToCells( TIteratorType GeometricalObjectsBegin, TIteratorType GeometricalObjectsEnd )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetCellBoundingBox](../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/GetCellBoundingBox.md)

# Called by

- [GeometricalObjectsBins](../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins.md)