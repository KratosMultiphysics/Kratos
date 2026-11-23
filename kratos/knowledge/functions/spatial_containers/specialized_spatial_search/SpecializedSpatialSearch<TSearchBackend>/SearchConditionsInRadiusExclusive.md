---
type: C++ Method
title: SearchConditionsInRadiusExclusive
resource: spatial_containers/specialized_spatial_search.cpp#L202-L269
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

`void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsInRadiusExclusive( const ConditionsContainerType& rStructureConditions, const ConditionsContainerType& rInputConditions, const RadiusArrayType& rRadius, VectorResultConditionsContainerType& rResults, VectorDistanceType& rResultsDistance )`

# Calls

- [PrepareSearch](../../../../functions/utilities/search_utilities/PrepareSearch.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [ParallelSearch](../../../../functions/utilities/search_utilities/ParallelSearch.md)