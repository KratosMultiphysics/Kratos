---
type: Python Function
title: simple_all_reduce_check
resource: mpi/tests/test_mpi_data_communicator_python.py#L315-L333
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/tests/test_mpi_data_communicator_python/check_value_in_rank
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/test_mpi_data_communicator_python/TestMPIDataCommunicatorPython/test_CastingTypes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def simple_all_reduce_check(method, ref_value: float):`

# Calls

- [check_value_in_rank](../../../../functions/mpi/tests/test_mpi_data_communicator_python/check_value_in_rank.md)

# Called by

- [test_CastingTypes](../../../../functions/mpi/tests/test_mpi_data_communicator_python/TestMPIDataCommunicatorPython/test_CastingTypes.md)