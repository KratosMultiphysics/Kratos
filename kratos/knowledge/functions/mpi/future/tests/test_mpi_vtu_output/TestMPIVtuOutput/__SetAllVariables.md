---
type: Python Method
title: __SetAllVariables
resource: mpi/future/tests/test_mpi_vtu_output.py#L154-L177
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__SetVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__InitializeModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __SetAllVariables(model_part: Kratos.ModelPart, offset = 0):`

# Calls

- [__SetVariable](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__SetVariable.md)

# Called by

- [__InitializeModelPart](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__InitializeModelPart.md)
- [__OutputTest](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest.md)