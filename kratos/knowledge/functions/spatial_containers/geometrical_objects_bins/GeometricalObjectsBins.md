---
type: C++ Function
title: GeometricalObjectsBins
resource: spatial_containers/geometrical_objects_bins.h#L75-L93
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GetTotalNumberOfCells
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/AddObjectsToCells
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GeometricalObjectsBins( TIteratorType GeometricalObjectsBegin, TIteratorType GeometricalObjectsEnd, const double Tolerance = 1e-12 )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetTotalNumberOfCells](../../../functions/spatial_containers/geometrical_objects_bins/GetTotalNumberOfCells.md)
- [AddObjectsToCells](../../../functions/spatial_containers/geometrical_objects_bins/AddObjectsToCells.md)