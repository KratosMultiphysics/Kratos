---
type: Python Function
title: check_value_in_rank
resource: mpi/tests/test_mpi_data_communicator_python.py#L284-L294
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/test_mpi_data_communicator_python/simple_reduce_check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_data_communicator_python/simple_all_reduce_check
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def check_value_in_rank(ref_value, value, ranks_to_check: 'list[int]'):`

# Calls

- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)

# Called by

- [simple_reduce_check](../../../../functions/mpi/tests/test_mpi_data_communicator_python/simple_reduce_check.md)
- [simple_all_reduce_check](../../../../functions/mpi/tests/test_mpi_data_communicator_python/simple_all_reduce_check.md)