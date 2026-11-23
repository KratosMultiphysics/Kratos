---
type: Python Method
title: __OutputTest
resource: mpi/future/tests/test_mpi_vtu_output.py#L184-L385
visibility: private
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
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__SetAllVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check_gauss
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/test_OutputASCII
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/test_OutputBinary
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __OutputTest(self, output_type: str):`

# Calls

- [VtuOutput](../../../../../../functions/input_output/vtu_output/VtuOutput/VtuOutput.md)
- [AddVariable](../../../../../../functions/input_output/vtu_output/VtuOutput/AddVariable.md)
- [AddIntegrationPointVariable](../../../../../../functions/input_output/vtu_output/VtuOutput/AddIntegrationPointVariable.md)
- [__SetAllVariables](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__SetAllVariables.md)
- [check](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/check.md)
- [check_gauss](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/check_gauss.md)
- [GetNumberOfNodes](../../../../../../functions/tests/test_vtu_output/TestVtuOutput/GetNumberOfNodes.md)

# Called by

- [test_OutputASCII](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/test_OutputASCII.md)
- [test_OutputBinary](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/test_OutputBinary.md)