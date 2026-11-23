---
type: C++ Method
title: VariableTensorAdaptor
resource: tensor_adaptors/variable_tensor_adaptor.cpp#L25-L40
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`VariableTensorAdaptor::VariableTensorAdaptor( TContainerPointerType pContainer, VariablePointerType pVariable) : mpVariable(pVariable)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)
- [DataShape](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape.md)