---
type: C++ Method
title: ComputeCommunicationPlanSubModelPart
resource: mpi/utilities/parallel_fill_communicator.cpp#L403-L522
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/SetLocalMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/SetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/InitializeParallelCommunicationMeshes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelFillCommunicator::ComputeCommunicationPlanSubModelPart( ModelPart& rSubModelPart, const std::vector<int>& rColors )`

# Calls

- [SetLocalMesh](../../../../../functions/sources/communicator/Communicator/SetLocalMesh.md)
- [SetNumberOfColors](../../../../../functions/sources/communicator/Communicator/SetNumberOfColors.md)

# Called by

- [InitializeParallelCommunicationMeshes](../../../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/InitializeParallelCommunicationMeshes.md)