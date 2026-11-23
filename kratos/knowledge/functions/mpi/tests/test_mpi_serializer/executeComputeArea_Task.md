---
type: Python Function
title: executeComputeArea_Task
resource: mpi/tests/test_mpi_serializer.py#L27-L58
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/sum
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/test_mpi_serializer/TestMPISerializer/testCalculateNodalArea
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def executeComputeArea_Task(pickled_model): # Unpickling model`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [sum](../../../../functions/includes/amatrix_interface/sum.md)

# Called by

- [testCalculateNodalArea](../../../../functions/mpi/tests/test_mpi_serializer/TestMPISerializer/testCalculateNodalArea.md)