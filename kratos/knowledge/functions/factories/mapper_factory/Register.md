---
type: C++ Function
title: Register
resource: factories/mapper_factory.h#L104-L112
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/factories/mapper_factory/GetApplicationPath
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry/Registry/GetCurrentSource
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void Register(const std::string& rMapperName, typename Mapper<TSparseSpace, TDenseSpace>::Pointer pMapperPrototype)`

# Calls

- [GetApplicationPath](../../../functions/factories/mapper_factory/GetApplicationPath.md)
- [GetCurrentSource](../../../functions/sources/registry/Registry/GetCurrentSource.md)