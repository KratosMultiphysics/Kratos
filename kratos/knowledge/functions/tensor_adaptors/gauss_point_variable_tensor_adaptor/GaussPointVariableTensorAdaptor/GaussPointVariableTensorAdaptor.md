---
type: C++ Method
title: GaussPointVariableTensorAdaptor
resource: tensor_adaptors/gauss_point_variable_tensor_adaptor.cpp#L29-L51
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

`GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor( TContainerPointerType pContainer, VariablePointerType pVariable, ProcessInfo::Pointer pProcessInfo) : mpVariable(pVariable), mpProcessInfo(pProcessInfo)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)
- [DataShape](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape.md)