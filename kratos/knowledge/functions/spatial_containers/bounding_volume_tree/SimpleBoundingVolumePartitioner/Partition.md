---
type: C++ Method
title: Partition
resource: spatial_containers/bounding_volume_tree.cpp#L68-L132
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/bounding_volume_tree/kDOP/GetLongestAxis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bounding_volume_tree/BoundingVolumePartitioner/ComputeCentroid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleBoundingVolumePartitioner::Partition(ConditionsContainerType& rAllConditions, const kDOP& rBoundingVolume, ConditionsContainerType& rOutputSet1, ConditionsContainerType& rOutputSet2)`

# Calls

- [GetLongestAxis](../../../../functions/spatial_containers/bounding_volume_tree/kDOP/GetLongestAxis.md)
- [ComputeCentroid](../../../../functions/spatial_containers/bounding_volume_tree/BoundingVolumePartitioner/ComputeCentroid.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)