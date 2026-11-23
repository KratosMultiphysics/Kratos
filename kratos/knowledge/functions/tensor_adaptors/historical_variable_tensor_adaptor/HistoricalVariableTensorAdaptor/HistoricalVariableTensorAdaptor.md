---
type: C++ Method
title: HistoricalVariableTensorAdaptor
resource: tensor_adaptors/historical_variable_tensor_adaptor.cpp#L57-L74
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

`HistoricalVariableTensorAdaptor::HistoricalVariableTensorAdaptor( ModelPart::NodesContainerType::Pointer pContainer, VariablePointerType pVariable, const int StepIndex) : mpVariable(pVariable), mStepIndex(StepIndex)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)
- [DataShape](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape.md)