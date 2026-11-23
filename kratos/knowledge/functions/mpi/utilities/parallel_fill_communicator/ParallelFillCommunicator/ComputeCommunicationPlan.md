---
type: C++ Method
title: ComputeCommunicationPlan
resource: mpi/utilities/parallel_fill_communicator.cpp#L174-L321
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Recv
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Send
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/InitializeParallelCommunicationMeshes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelFillCommunicator::ComputeCommunicationPlan(ModelPart& rModelPart)`

# Calls

- [SolutionStepsDataHas](../../../../../functions/includes/node/SolutionStepsDataHas.md)
- [Recv](../../../../../functions/includes/data_communicator/Recv.md)
- [Send](../../../../../functions/includes/data_communicator/Send.md)
- [Broadcast](../../../../../functions/includes/data_communicator/Broadcast.md)
- [InitializeParallelCommunicationMeshes](../../../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/InitializeParallelCommunicationMeshes.md)