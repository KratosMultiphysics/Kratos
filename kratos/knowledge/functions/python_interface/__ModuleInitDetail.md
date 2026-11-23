---
type: Python Function
title: __ModuleInitDetail
resource: python_interface/__init__.py#L70-L126
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/IsMPIAvailable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python/kratos_mpi_python/InitializeMPIParallelRun
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __ModuleInitDetail():`

# Calls

- [IsMPIAvailable](../../functions/python_scripts/kratos_utilities/IsMPIAvailable.md)
- [InitializeMPIParallelRun](../../functions/mpi/python/kratos_mpi_python/InitializeMPIParallelRun.md)