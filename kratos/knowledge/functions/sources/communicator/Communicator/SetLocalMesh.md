---
type: C++ Method
title: SetLocalMesh
resource: sources/communicator.cpp#L173-L176
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/DuplicateCommunicatorData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlanSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/ModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Communicator::SetLocalMesh(MeshType::Pointer pGivenMesh)`

# Called by

- [DuplicateCommunicatorData](../../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/DuplicateCommunicatorData.md)
- [ComputeCommunicationPlanSubModelPart](../../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlanSubModelPart.md)
- [ModelPart](../../../../functions/sources/model_part/ModelPart/ModelPart.md)