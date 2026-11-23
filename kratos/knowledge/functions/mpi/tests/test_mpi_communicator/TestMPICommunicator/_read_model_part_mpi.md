---
type: Python Method
title: _read_model_part_mpi
resource: mpi/tests/test_mpi_communicator.py#L16-L30
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_assemble_variable_in_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_assemble_variable_in_sub_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorRankSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorReduction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _read_model_part_mpi(self,main_model_part):`

# Calls

- [HasSubModelPart](../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [test_assemble_variable_in_model_part](../../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_assemble_variable_in_model_part.md)
- [test_assemble_variable_in_sub_model_part](../../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_assemble_variable_in_sub_model_part.md)
- [testCommunicatorRankSize](../../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorRankSize.md)
- [testCommunicatorReduction](../../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorReduction.md)
- [test_GlobalNumberOf_Methods](../../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods.md)