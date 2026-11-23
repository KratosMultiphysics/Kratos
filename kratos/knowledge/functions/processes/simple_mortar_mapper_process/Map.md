---
type: C++ Function
title: Map
resource: processes/simple_mortar_mapper_process.h#L257-L269
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python/add_mapper_to_python/MapWithoutOptionsScalar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mapper_to_python/MapWithoutOptionsVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mapper_to_python/MapWithOptionsScalar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mapper_to_python/MapWithOptionsVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Map( const TVarType& rOriginVariable, const TVarType& rDestinationVariable, const Flags Flag = Flags() )`

# Called by

- [MapWithoutOptionsScalar](../../../functions/python/add_mapper_to_python/MapWithoutOptionsScalar.md)
- [MapWithoutOptionsVector](../../../functions/python/add_mapper_to_python/MapWithoutOptionsVector.md)
- [MapWithOptionsScalar](../../../functions/python/add_mapper_to_python/MapWithOptionsScalar.md)
- [MapWithOptionsVector](../../../functions/python/add_mapper_to_python/MapWithOptionsVector.md)