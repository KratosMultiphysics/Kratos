---
type: Python Method
title: GetUnstructuredGridName
resource: tests/test_vtu_output.py#L436-L448
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetUnstructuredGridName(unstructured_grid: 'tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]', prefix: str, step: int, data_communicator: Kratos.DataCommunicator, entity_suffix = "s") -> str:`

# Calls

- [FullName](../../../../functions/includes/model_part/FullName.md)

# Called by

- [test_PointVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)