---
type: Python Method
title: CheckGaussVtuFile
resource: tests/test_vtu_output.py#L582-L612
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check_gauss
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

`def CheckGaussVtuFile(test_class: kratos_unittest.TestCase, model_part: Kratos.ModelPart, container, prefix: str, step_id: int, output_type: str, data_communicator: Kratos.DataCommunicator, point_fields):`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [CheckVtuFile](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile.md)
- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)

# Called by

- [check_gauss](../../../../functions/mpi/future/tests/test_mpi_vtu_output/check_gauss.md)
- [test_PointVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)