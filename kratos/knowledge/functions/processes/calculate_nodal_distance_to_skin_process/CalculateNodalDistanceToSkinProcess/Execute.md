---
type: C++ Method
title: Execute
resource: processes/calculate_nodal_distance_to_skin_process.cpp#L92-L154
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/spatial_search_result/SpatialSearchResult/GetDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetNonHistoricalVariableToZero
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateNodalDistanceToSkinProcess::Execute()`

# Calls

- [GetDistance](../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/GetDistance.md)
- [SetNonHistoricalVariableToZero](../../../../functions/utilities/variable_utils/SetNonHistoricalVariableToZero.md)