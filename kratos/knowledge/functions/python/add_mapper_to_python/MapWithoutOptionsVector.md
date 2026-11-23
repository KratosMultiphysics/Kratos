---
type: C++ Function
title: MapWithoutOptionsVector
resource: python/add_mapper_to_python.h#L66-L72
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/simple_mortar_mapper_process/Map
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void MapWithoutOptionsVector(Mapper<TSparseSpace, TDenseSpace>& dummy, const Variable< array_1d<double, 3> >& origin_variable, const Variable< array_1d<double, 3> >& destination_variable)`

# Calls

- [Map](../../../functions/processes/simple_mortar_mapper_process/Map.md)