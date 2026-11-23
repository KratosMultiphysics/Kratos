---
type: C++ Function
title: DeepCopyEntities
resource: utilities/auxiliar_model_part_utilities.h#L684-L711
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DeepCopyEntities( ModelPart& rModelPart, TClassContainer& rEntities, TReferenceClassContainer& rReferenceEntities, std::unordered_map<Geometry<Node>::Pointer,Geometry<Node>::Pointer>& rGeometryPointerDatabase )`

# Calls

- [GetData](../../../functions/sources/gid_io/GidIOBase/GetData.md)