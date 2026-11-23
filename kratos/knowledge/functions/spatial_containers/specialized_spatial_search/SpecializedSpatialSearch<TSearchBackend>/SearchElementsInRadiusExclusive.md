---
type: C++ Method
title: SearchElementsInRadiusExclusive
resource: spatial_containers/specialized_spatial_search.cpp#L26-L93
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/search_utilities/PrepareSearch
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/ParallelSearch
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SpecializedSpatialSearch<TSearchBackend>::SearchElementsInRadiusExclusive( const ElementsContainerType& rStructureElements, const ElementsContainerType& rInputElements, const RadiusArrayType& rRadius, VectorResultElementsContainerType& rResults, VectorDistanceType& rResultsDistance )`

# Calls

- [PrepareSearch](../../../../functions/utilities/search_utilities/PrepareSearch.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [ParallelSearch](../../../../functions/utilities/search_utilities/ParallelSearch.md)