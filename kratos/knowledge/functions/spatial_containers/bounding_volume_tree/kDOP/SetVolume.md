---
type: C++ Method
title: SetVolume
resource: spatial_containers/bounding_volume_tree.h#L152-L162
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/bounding_volume_tree/kDOP/GetType
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spatial_containers/bounding_volume_tree/BoundingVolumeTree/UpdateTree
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetVolume(const kDOP& rBV1, const kDOP& rBV2)`

# Calls

- [GetType](../../../../functions/spatial_containers/bounding_volume_tree/kDOP/GetType.md)

# Called by

- [UpdateTree](../../../../functions/spatial_containers/bounding_volume_tree/BoundingVolumeTree/UpdateTree.md)