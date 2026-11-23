---
type: C++ Function
title: CreateMapper
resource: factories/mapper_factory.h#L65-L102
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry_item/RegistryItem/KeyConstBegin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry_item/RegistryItem/KeyConstEnd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename Mapper<TSparseSpace, TDenseSpace>::Pointer CreateMapper( ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination, Parameters MapperSettings)`

# Calls

- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RemoveValue](../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [KeyConstBegin](../../../functions/sources/registry_item/RegistryItem/KeyConstBegin.md)
- [KeyConstEnd](../../../functions/sources/registry_item/RegistryItem/KeyConstEnd.md)