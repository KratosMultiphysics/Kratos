---
type: C++ Method
title: pGetObject
resource: spatial_containers/point_object.cpp#L70-L73
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/FillDatabase
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/ParallelSearch
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename TObject::Pointer PointObject<TObject>::pGetObject() const`

# Called by

- [CalculateDistance](../../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistance.md)
- [FillDatabase](../../../../functions/processes/simple_mortar_mapper_process/FillDatabase.md)
- [ParallelSearch](../../../../functions/utilities/search_utilities/ParallelSearch.md)