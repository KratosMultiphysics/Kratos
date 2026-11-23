---
type: C++ Function
title: PointerDistance
resource: spatial_containers/search_structure.h#L91-L95
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/spatial_containers/bins_dynamic/BinsDynamic/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bins_dynamic/BinsDynamic/Construct
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bins_dynamic_objects/BinsObjectDynamic/BinsObjectDynamic
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bins_static/Bins/GenerateBins
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bins_static/Bins/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bins_static/Bins/Construct
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bins_static_objects/BinsObjectStatic/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bins_static_objects/BinsObjectStatic/CalculateBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bucket/Bucket/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/kd_tree/KDTreePartition/Partition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/kd_tree/KDTreePartition/Construct
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/MaxSpread
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/Construct
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/kd_tree/KDTreePartitionMidPointSplit/Construct
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree/OCTreePartition/OCTreePartition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree/OCTreePartition/Construct
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/tree/Tree/Tree
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t PointerDistance( TPointerType const& PointerBegin, TPointerType const& PointerEnd )`

# Called by

- [PrintData](../../../functions/spatial_containers/bins_dynamic/BinsDynamic/PrintData.md)
- [Construct](../../../functions/spatial_containers/bins_dynamic/BinsDynamic/Construct.md)
- [BinsObjectDynamic](../../../functions/spatial_containers/bins_dynamic_objects/BinsObjectDynamic/BinsObjectDynamic.md)
- [GenerateBins](../../../functions/spatial_containers/bins_static/Bins/GenerateBins.md)
- [PrintData](../../../functions/spatial_containers/bins_static/Bins/PrintData.md)
- [Construct](../../../functions/spatial_containers/bins_static/Bins/Construct.md)
- [PrintData](../../../functions/spatial_containers/bins_static_objects/BinsObjectStatic/PrintData.md)
- [CalculateBoundingBox](../../../functions/spatial_containers/bins_static_objects/BinsObjectStatic/CalculateBoundingBox.md)
- [PrintData](../../../functions/spatial_containers/bucket/Bucket/PrintData.md)
- [Partition](../../../functions/spatial_containers/kd_tree/KDTreePartition/Partition.md)
- [Construct](../../../functions/spatial_containers/kd_tree/KDTreePartition/Construct.md)
- [MaxSpread](../../../functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/MaxSpread.md)
- [Construct](../../../functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/Construct.md)
- [Construct](../../../functions/spatial_containers/kd_tree/KDTreePartitionMidPointSplit/Construct.md)
- [OCTreePartition](../../../functions/spatial_containers/octree/OCTreePartition/OCTreePartition.md)
- [Construct](../../../functions/spatial_containers/octree/OCTreePartition/Construct.md)
- [Tree](../../../functions/spatial_containers/tree/Tree/Tree.md)