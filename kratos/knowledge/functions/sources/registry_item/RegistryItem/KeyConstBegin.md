---
type: C++ Method
title: KeyConstBegin
resource: sources/registry_item.cpp#L63-L66
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/factories/mapper_factory/CreateMapper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/factories/mapper_factory/GetRegisteredMapperNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_registry_to_python/AddRegistryToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_application/KratosApplication/DeregisterVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_application/KratosApplication/DeregisterMappers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_registry/CheckKeyIteration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`RegistryItem::KeyReturnConstIterator RegistryItem::KeyConstBegin() const`

# Called by

- [CreateMapper](../../../../functions/factories/mapper_factory/CreateMapper.md)
- [GetRegisteredMapperNames](../../../../functions/factories/mapper_factory/GetRegisteredMapperNames.md)
- [AddRegistryToPython](../../../../functions/python/add_registry_to_python/AddRegistryToPython.md)
- [DeregisterVariables](../../../../functions/sources/kratos_application/KratosApplication/DeregisterVariables.md)
- [DeregisterMappers](../../../../functions/sources/kratos_application/KratosApplication/DeregisterMappers.md)
- [CheckKeyIteration](../../../../functions/tests/cpp_tests/sources/test_registry/CheckKeyIteration.md)