---
type: C++ Function
title: ParallelSearch
resource: utilities/search_utilities.h#L551-L587
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/point_object/PointObject<TObject>/pGetObject
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchElementsInRadiusExclusive
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchNodesInRadiusExclusive
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchConditionsInRadiusExclusive
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ParallelSearch( const TContainer& rInput, const RadiusArrayType& rRadius, TSpatialContainer& rSearch, TResultType& rResults, VectorDistanceType& rResultsDistance, const int AllocationSize = 1000 )`

# Calls

- [pGetObject](../../../functions/spatial_containers/point_object/PointObject<TObject>/pGetObject.md)

# Called by

- [SearchElementsInRadiusExclusive](../../../functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchElementsInRadiusExclusive.md)
- [SearchNodesInRadiusExclusive](../../../functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchNodesInRadiusExclusive.md)
- [SearchConditionsInRadiusExclusive](../../../functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchConditionsInRadiusExclusive.md)