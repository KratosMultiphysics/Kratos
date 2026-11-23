---
type: Python Module
title: distributed_restart_utility
resource: mpi/python_scripts/distributed_restart_utility.py#L1-L65
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-mpi-as-kratosmpi
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-restart-utility
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [DistributedRestartUtility](../../../classes/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility.md)
- [__init__](../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/__init__.md)
- [_GetFileLabelLoad](../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_GetFileLabelLoad.md)
- [_GetFileLabelSave](../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_GetFileLabelSave.md)
- [_ExecuteAfterLoad](../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_ExecuteAfterLoad.md)
- [_GetFileNamePattern](../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_GetFileNamePattern.md)
- [_UpdateRestartFilesMap](../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_UpdateRestartFilesMap.md)
- [_GetSerializerFlags](../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_GetSerializerFlags.md)
- [_GetDefaultParameters](../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_GetDefaultParameters.md)

# Imports

- `KratosMultiphysics`
- `KratosMultiphysics.mpi as KratosMPI`
- `KratosMultiphysics.restart_utility`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)