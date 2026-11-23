---
type: Python Module
title: distributed_import_model_part_utility
resource: mpi/python_scripts/distributed_import_model_part_utility.py#L1-L144
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
  - target: external/pathlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/time
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-metisapplication-as-kratosmetis
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-mpi-distributed-restart-utility
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [DistributedImportModelPartUtility](../../../classes/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility.md)
- [__init__](../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/__init__.md)
- [ExecutePartitioningAndReading](../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ExecutePartitioningAndReading.md)
- [ImportModelPart](../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart.md)
- [CreateCommunicators](../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/CreateCommunicators.md)

# Imports

- `KratosMultiphysics`
- `KratosMultiphysics.mpi as KratosMPI`
- `pathlib`
- `time`
- `KratosMultiphysics.MetisApplication as KratosMetis`
- `KratosMultiphysics.mpi.distributed_restart_utility`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)