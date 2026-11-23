---
type: Python Method
title: GetNumberOfNodes
resource: tests/test_vtu_output.py#L461-L466
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetNumberOfNodes(container: 'typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]') -> int:`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [__OutputTest](../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest.md)
- [test_PointVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)