---
type: Python Method
title: __TestFlagsTensorAdaptor
resource: tests/test_tensor_adaptors.py#L847-L881
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodalFlagsTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ConditionFlagsTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ElementFlagsTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __TestFlagsTensorAdaptor(self, container):`

# Calls

- [IsDefined](../../../../functions/containers/flags/IsDefined.md)

# Called by

- [test_NodalFlagsTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodalFlagsTensorAdaptor.md)
- [test_ConditionFlagsTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ConditionFlagsTensorAdaptor.md)
- [test_ElementFlagsTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ElementFlagsTensorAdaptor.md)