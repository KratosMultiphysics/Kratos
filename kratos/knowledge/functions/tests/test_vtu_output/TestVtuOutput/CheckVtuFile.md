---
type: Python Method
title: CheckVtuFile
resource: tests/test_vtu_output.py#L477-L537
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckDataArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CheckVtuFile( test_class: kratos_unittest.TestCase, vtu_file_name: str, number_of_points: int, number_of_cells: 'typing.Optional[int]', data_format : str, point_data_fields: 'dict[str, tuple[int, str]]', cell_data_fields: 'dict[str, tuple[int, str]]'):`

# Calls

- [CheckDataArray](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckDataArray.md)
- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)

# Called by

- [check](../../../../functions/mpi/future/tests/test_mpi_vtu_output/check.md)
- [test_PointVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)
- [CheckGaussVtuFile](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile.md)