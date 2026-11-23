---
type: C++ Method
title: InitializeParallelCommunicationMeshes
resource: mpi/utilities/parallel_fill_communicator.cpp#L323-L401
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/SetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/GenerateMeshes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlanSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlan
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelFillCommunicator::InitializeParallelCommunicationMeshes( ModelPart& rModelPart, const std::vector<int>& rColors, const int MyRank )`

# Calls

- [SetNumberOfColors](../../../../../functions/sources/communicator/Communicator/SetNumberOfColors.md)
- [GenerateMeshes](../../../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/GenerateMeshes.md)
- [ComputeCommunicationPlanSubModelPart](../../../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlanSubModelPart.md)

# Called by

- [ComputeCommunicationPlan](../../../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlan.md)