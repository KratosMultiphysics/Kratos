---
type: C++ Method
title: InsertPoint
resource: spatial_containers/bounding_volume_tree.h#L119-L129
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/bounding_volume_tree/kDOP/NormalCoordinate
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spatial_containers/bounding_volume_tree/InsertGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InsertPoint(const double& rX, const double& rY, const double& rZ)`

# Calls

- [NormalCoordinate](../../../../functions/spatial_containers/bounding_volume_tree/kDOP/NormalCoordinate.md)

# Called by

- [InsertGeometry](../../../../functions/spatial_containers/bounding_volume_tree/InsertGeometry.md)