---
type: Python Function
title: check_pvtu
resource: mpi/future/tests/test_mpi_vtu_output.py#L273-L324
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
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
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check_gauss
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def check_pvtu(model_part: Kratos.ModelPart, suffix: str, step_id: int, point_fields, cell_fields):`

# Calls

- [FullName](../../../../../functions/includes/model_part/FullName.md)
- [CheckDataArray](../../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckDataArray.md)
- [DeleteFileIfExisting](../../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)

# Called by

- [check](../../../../../functions/mpi/future/tests/test_mpi_vtu_output/check.md)
- [check_gauss](../../../../../functions/mpi/future/tests/test_mpi_vtu_output/check_gauss.md)