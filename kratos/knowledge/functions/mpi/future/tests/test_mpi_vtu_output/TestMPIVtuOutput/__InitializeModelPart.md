---
type: Python Method
title: __InitializeModelPart
resource: mpi/future/tests/test_mpi_vtu_output.py#L73-L130
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__CreateEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__SetAllVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/setUpClass
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __InitializeModelPart(model_part: Kratos.ModelPart): # Add variables.`

# Calls

- [__CreateEntities](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__CreateEntities.md)
- [__SetAllVariables](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__SetAllVariables.md)
- [GetGeometry](../../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [setUpClass](../../../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/setUpClass.md)