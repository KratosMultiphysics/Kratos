---
type: Python Function
title: CheckIfApplicationsAvailable
resource: python_scripts/kratos_utilities.py#L101-L105
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/GetListOfAvailableApplications
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/python_scripts/distributed_gid_output_process/DistributedGiDOutputProcess/_CreateCuttingUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/eigen_solver_factory/ConstructSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_linear_solver_factory/CreateFastestAvailableDirectLinearSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CheckIfApplicationsAvailable(*application_names):`

# Calls

- [GetListOfAvailableApplications](../../../functions/python_scripts/kratos_utilities/GetListOfAvailableApplications.md)

# Called by

- [_CreateCuttingUtility](../../../functions/mpi/python_scripts/distributed_gid_output_process/DistributedGiDOutputProcess/_CreateCuttingUtility.md)
- [ConstructSolver](../../../functions/python_scripts/eigen_solver_factory/ConstructSolver.md)
- [CreateFastestAvailableDirectLinearSolver](../../../functions/python_scripts/python_linear_solver_factory/CreateFastestAvailableDirectLinearSolver.md)