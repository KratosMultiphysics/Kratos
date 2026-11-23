---
type: Python Method
title: __GetUnstructuredGridList
resource: tests/test_vtu_output.py#L615-L632
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __GetUnstructuredGridList(unstructured_grid_list: 'list[tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]]', model_part: Kratos.ModelPart, data_communicator: Kratos.DataCommunicator, recursively: bool) -> None:`

# Called by

- [GetUnstructuredGridList](../../../../functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridList.md)