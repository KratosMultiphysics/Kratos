---
type: Python Function
title: check_gauss
resource: mpi/future/tests/test_mpi_vtu_output.py#L336-L356
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile
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

`def check_gauss(model_part: Kratos.ModelPart, container, step_id: int, point_fields):`

# Calls

- [CheckGaussVtuFile](../../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile.md)
- [DeleteFileIfExisting](../../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)
- [check_pvtu](../../../../../functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu.md)

# Called by

- [__OutputTest](../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest.md)