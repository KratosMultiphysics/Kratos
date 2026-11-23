---
type: C++ Method
title: FixityTensorAdaptor
resource: tensor_adaptors/fixity_tensor_adaptor.cpp#L38-L57
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

`FixityTensorAdaptor::FixityTensorAdaptor( ContainerPointerType pContainer, const std::vector<const Variable<double>*>& rDofsVarPointerList) : TensorAdaptor<bool>() , mDofsVarPointerList(rDofsVarPointerList)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)