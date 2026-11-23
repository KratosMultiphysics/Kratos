---
type: Python Function
title: DeleteDirectoryIfExistingAndEmpty
resource: python_scripts/kratos_utilities.py#L26-L32
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/tearDownClass
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput2D/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput3D/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/tearDownClass
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def DeleteDirectoryIfExistingAndEmpty(directory_name):`

# Calls

- [DeleteDirectoryIfExisting](../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)

# Called by

- [tearDownClass](../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/tearDownClass.md)
- [tearDown](../../../functions/tests/test_vtu_output/TestVtuOutput2D/tearDown.md)
- [tearDown](../../../functions/tests/test_vtu_output/TestVtuOutput3D/tearDown.md)
- [tearDownClass](../../../functions/tests/test_vtu_output/TestVtuOutput/tearDownClass.md)