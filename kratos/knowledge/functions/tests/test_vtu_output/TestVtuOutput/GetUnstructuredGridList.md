---
type: Python Method
title: GetUnstructuredGridList
resource: tests/test_vtu_output.py#L430-L433
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_vtu_output/TestVtuOutput/__GetUnstructuredGridList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_GetOutputContainerList
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

`def GetUnstructuredGridList(model_part: Kratos.ModelPart, data_communicator: Kratos.DataCommunicator, recursively: bool) -> 'list[tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]]':`

# Calls

- [__GetUnstructuredGridList](../../../../functions/tests/test_vtu_output/TestVtuOutput/__GetUnstructuredGridList.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)

# Called by

- [test_GetOutputContainerList](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_GetOutputContainerList.md)
- [test_PointVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)