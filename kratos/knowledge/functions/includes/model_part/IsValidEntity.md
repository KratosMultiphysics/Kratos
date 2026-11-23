---
type: C++ Function
title: IsValidEntity
resource: includes/model_part.h#L1993-L2049
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part/ModelPart/AddGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static bool IsValidEntity(const TContainerType& r_entities, const typename TContainerType::value_type * r_entity, std::string_view rModelPartName, std::string_view rRootModelPartName)`

# Calls

- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [AddGeometry](../../../functions/sources/model_part/ModelPart/AddGeometry.md)