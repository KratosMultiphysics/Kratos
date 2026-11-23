---
type: Python Method
title: AddItem
resource: python_interface/python_registry.py#L218-L229
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/python_registry/PythonRegistry/__AddItemToAll
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry/PythonRegistry/__AddItem
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_interface/python_registry/register_wrapper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry_utilities/__RegisterItemList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testHasItemPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testNumberOfItems
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testAddItem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testAddItemRepeatedInAllBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testAddItemAlreadyRegisteredInCpp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testAddItemAlreadyRegisteredInPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testRemoveItemPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testGetItemPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_python_registry/TestPythonRegistry/testKeys
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def AddItem(self, Name, ModuleName): # Add current item # First check if it is already registered in the c++ and Python registries # If not registered, add it to the Python registry (note that we do an "All" registry as well as the module one)`

# Calls

- [__AddItemToAll](../../../../functions/python_interface/python_registry/PythonRegistry/__AddItemToAll.md)
- [__AddItem](../../../../functions/python_interface/python_registry/PythonRegistry/__AddItem.md)

# Called by

- [register_wrapper](../../../../functions/python_interface/python_registry/register_wrapper.md)
- [__RegisterItemList](../../../../functions/python_interface/python_registry_utilities/__RegisterItemList.md)
- [testHasItemPython](../../../../functions/tests/test_python_registry/TestPythonRegistry/testHasItemPython.md)
- [testNumberOfItems](../../../../functions/tests/test_python_registry/TestPythonRegistry/testNumberOfItems.md)
- [testAddItem](../../../../functions/tests/test_python_registry/TestPythonRegistry/testAddItem.md)
- [testAddItemRepeatedInAllBlock](../../../../functions/tests/test_python_registry/TestPythonRegistry/testAddItemRepeatedInAllBlock.md)
- [testAddItemAlreadyRegisteredInCpp](../../../../functions/tests/test_python_registry/TestPythonRegistry/testAddItemAlreadyRegisteredInCpp.md)
- [testAddItemAlreadyRegisteredInPython](../../../../functions/tests/test_python_registry/TestPythonRegistry/testAddItemAlreadyRegisteredInPython.md)
- [testRemoveItemPython](../../../../functions/tests/test_python_registry/TestPythonRegistry/testRemoveItemPython.md)
- [testGetItemPython](../../../../functions/tests/test_python_registry/TestPythonRegistry/testGetItemPython.md)
- [testKeys](../../../../functions/tests/test_python_registry/TestPythonRegistry/testKeys.md)