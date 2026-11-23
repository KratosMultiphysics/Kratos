---
type: Python Method
title: test_CellVariableAddition
resource: tests/test_vtu_output.py#L375-L427
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/VtuOutput/VtuOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/AddVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/AddIntegrationPointVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckVtmFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckSeriesFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_CellVariableAddition(self):`

# Calls

- [VtuOutput](../../../../functions/input_output/vtu_output/VtuOutput/VtuOutput.md)
- [AddVariable](../../../../functions/input_output/vtu_output/VtuOutput/AddVariable.md)
- [AddIntegrationPointVariable](../../../../functions/input_output/vtu_output/VtuOutput/AddIntegrationPointVariable.md)
- [GetUnstructuredGridList](../../../../functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridList.md)
- [GetUnstructuredGridName](../../../../functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridName.md)
- [CheckVtuFile](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile.md)
- [GetNumberOfNodes](../../../../functions/tests/test_vtu_output/TestVtuOutput/GetNumberOfNodes.md)
- [CheckGaussVtuFile](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile.md)
- [CheckVtmFile](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtmFile.md)
- [CheckSeriesFile](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckSeriesFile.md)