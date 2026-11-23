---
type: Python Method
title: test_GlobalNumberOf_Methods
resource: mpi/tests/test_mpi_communicator.py#L181-L219
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/_read_model_part_mpi
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/VariableUtils/AddDofsList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_GlobalNumberOf_Methods(self):`

# Calls

- [_read_model_part_mpi](../../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/_read_model_part_mpi.md)
- [AddDofsList](../../../../../functions/utilities/variable_utils/VariableUtils/AddDofsList.md)
- [GlobalNumberOfNodes](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfNodes.md)
- [GlobalNumberOfElements](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)
- [GlobalNumberOfMasterSlaveConstraints](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfMasterSlaveConstraints.md)