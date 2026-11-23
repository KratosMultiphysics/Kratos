---
type: Python Method
title: tearDown
resource: mpi/tests/test_distributed_import_model_part_utility.py#L15-L25
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def tearDown(self): # Remove the Metis partitioning files # this can only be done after all processes arrived here! # In these tests not all ranks participate in the test (=> SubComm), # hence they could remove the files for the other ranks before they could read!`

# Calls

- [DeleteDirectoryIfExisting](../../../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)