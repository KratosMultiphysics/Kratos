---
type: Python Method
title: __TestVariableTensorAdaptorPartial
resource: tests/test_tensor_adaptors.py#L916-L1020
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodalVariableTensorAdaptorPartial
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ConditionVariableTensorAdaptorPartial
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ElementVariableTensorAdaptorPartial
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_PropertiesVariableTensorAdaptorPartial
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __TestVariableTensorAdaptorPartial(self, container):`

# Calls

- [DataShape](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape.md)

# Called by

- [test_NodalVariableTensorAdaptorPartial](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodalVariableTensorAdaptorPartial.md)
- [test_ConditionVariableTensorAdaptorPartial](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ConditionVariableTensorAdaptorPartial.md)
- [test_ElementVariableTensorAdaptorPartial](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ElementVariableTensorAdaptorPartial.md)
- [test_PropertiesVariableTensorAdaptorPartial](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_PropertiesVariableTensorAdaptorPartial.md)