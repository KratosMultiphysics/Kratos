---
type: C++ Method
title: FlagsTensorAdaptor
resource: tensor_adaptors/flags_tensor_adaptor.cpp#L34-L40
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FlagsTensorAdaptor::FlagsTensorAdaptor( TContainerPointerType pContainer, const Flags& rFlags) : BaseType(pContainer, Kratos::make_shared<Storage>(DenseVector<unsigned int>(1, pContainer->size()))), mFlags(rFlags)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)