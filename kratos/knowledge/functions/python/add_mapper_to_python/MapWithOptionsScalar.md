---
type: C++ Function
title: MapWithOptionsScalar
resource: python/add_mapper_to_python.h#L94-L100
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/simple_mortar_mapper_process/Map
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void MapWithOptionsScalar(Mapper<TSparseSpace, TDenseSpace>& dummy, const Variable<double>& origin_variable, const Variable<double>& destination_variable, Kratos::Flags MappingOptions)`

# Calls

- [Map](../../../functions/processes/simple_mortar_mapper_process/Map.md)