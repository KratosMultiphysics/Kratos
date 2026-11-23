---
type: C++ Function
title: GenerateSpatialSearchSolutionTuple
resource: python/add_search_strategies_to_python.cpp#L120-L129
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/add_search_strategies_to_python/GenerateListFromVectorOfVector
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_search_strategies_to_python/AddSearchStrategiesToPython
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::tuple<pybind11::list, pybind11::list> GenerateSpatialSearchSolutionTuple( const TClass& rResults, SpatialSearch::VectorDistanceType& rDistances )`

# Calls

- [GenerateListFromVectorOfVector](../../../functions/python/add_search_strategies_to_python/GenerateListFromVectorOfVector.md)

# Called by

- [AddSearchStrategiesToPython](../../../functions/python/add_search_strategies_to_python/AddSearchStrategiesToPython.md)