---
type: C++ Function
title: MapWithOptionsVector
resource: python/add_mapper_to_python.h#L103-L109
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/simple_mortar_mapper_process/Map
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void MapWithOptionsVector(Mapper<TSparseSpace, TDenseSpace>& dummy, const Variable< array_1d<double, 3> >& origin_variable, const Variable< array_1d<double, 3> >& destination_variable, Kratos::Flags MappingOptions)`

# Calls

- [Map](../../../functions/processes/simple_mortar_mapper_process/Map.md)