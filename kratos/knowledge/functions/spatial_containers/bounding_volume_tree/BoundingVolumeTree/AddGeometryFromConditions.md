---
type: C++ Method
title: AddGeometryFromConditions
resource: spatial_containers/bounding_volume_tree.h#L533-L540
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spatial_containers/bounding_volume_tree/BoundingVolumeTree/BuildTreeTopDown
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddGeometryFromConditions(ConditionsContainerType& rAllConditions)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [BuildTreeTopDown](../../../../functions/spatial_containers/bounding_volume_tree/BoundingVolumeTree/BuildTreeTopDown.md)