---
type: Python Function
title: check
resource: mpi/future/tests/test_mpi_vtu_output.py#L326-L334
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def check(model_part: Kratos.ModelPart, number_of_nodes: int, number_of_cells: int, step_id: int, suffix: str, point_fields, cell_fields):`

# Calls

- [FullName](../../../../../functions/includes/model_part/FullName.md)
- [CheckVtuFile](../../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile.md)
- [DeleteFileIfExisting](../../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)
- [check_pvtu](../../../../../functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu.md)

# Called by

- [__OutputTest](../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest.md)