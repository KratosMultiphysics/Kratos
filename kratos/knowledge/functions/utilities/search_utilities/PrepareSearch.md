---
type: C++ Function
title: PrepareSearch
resource: utilities/search_utilities.h#L449-L461
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/search_utilities/PrepareOutputSearch
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/PreparePointsSearch
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

`static std::vector<typename PointObject<typename TContainer::value_type>::Pointer> PrepareSearch( const TContainer& rStructure, const TContainer& rInput, TResultType& rResults, VectorDistanceType& rResultsDistance )`

# Calls

- [PrepareOutputSearch](../../../functions/utilities/search_utilities/PrepareOutputSearch.md)
- [PreparePointsSearch](../../../functions/utilities/search_utilities/PreparePointsSearch.md)

# Called by

- [SearchElementsInRadiusExclusive](../../../functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchElementsInRadiusExclusive.md)
- [SearchNodesInRadiusExclusive](../../../functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchNodesInRadiusExclusive.md)
- [SearchConditionsInRadiusExclusive](../../../functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch<TSearchBackend>/SearchConditionsInRadiusExclusive.md)