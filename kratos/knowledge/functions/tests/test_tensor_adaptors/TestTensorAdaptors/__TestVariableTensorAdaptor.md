---
type: Python Method
title: __TestVariableTensorAdaptor
resource: tests/test_tensor_adaptors.py#L895-L914
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__CheckValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodeVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ConditionVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ElementVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_PropertyVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __TestVariableTensorAdaptor(self, container):`

# Calls

- [DataShape](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape.md)
- [__CheckValues](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__CheckValues.md)

# Called by

- [test_NodeVariableTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodeVariableTensorAdaptor.md)
- [test_ConditionVariableTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ConditionVariableTensorAdaptor.md)
- [test_ElementVariableTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ElementVariableTensorAdaptor.md)
- [test_PropertyVariableTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_PropertyVariableTensorAdaptor.md)