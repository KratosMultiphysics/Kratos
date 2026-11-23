---
type: Python Method
title: LoadRestart
resource: python_scripts/restart_utility.py#L88-L117
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/restart_utility/RestartUtility/__GetFileNameLoad
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/exists
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_solver/PythonSolver/_ImportModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/__execute_restart_utility_load
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def LoadRestart(self, restart_file_name=""):`

# Calls

- [__GetFileNameLoad](../../../../functions/python_scripts/restart_utility/RestartUtility/__GetFileNameLoad.md)
- [exists](../../../../functions/sources/kratos_filesystem/exists.md)

# Called by

- [ImportModelPart](../../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart.md)
- [_ImportModelPart](../../../../functions/python_scripts/python_solver/PythonSolver/_ImportModelPart.md)
- [__execute_restart_utility_load](../../../../functions/tests/test_restart/TestRestart/__execute_restart_utility_load.md)