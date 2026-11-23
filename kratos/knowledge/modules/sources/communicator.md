---
type: C++ Module
title: communicator
resource: sources/communicator.cpp#L1-L672
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-parallel-environment-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [Communicator](../../functions/sources/communicator/Communicator/Communicator.md)
- [Communicator](../../functions/sources/communicator/Communicator/Communicator-2.md)
- [Communicator](../../functions/sources/communicator/Communicator/Communicator-3.md)
- [Create](../../functions/sources/communicator/Communicator/Create.md)
- [Create](../../functions/sources/communicator/Communicator/Create-2.md)
- [IsDistributed](../../functions/sources/communicator/Communicator/IsDistributed.md)
- [MyPID](../../functions/sources/communicator/Communicator/MyPID.md)
- [TotalProcesses](../../functions/sources/communicator/Communicator/TotalProcesses.md)
- [GlobalNumberOfNodes](../../functions/sources/communicator/Communicator/GlobalNumberOfNodes.md)
- [GlobalNumberOfElements](../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)
- [GlobalNumberOfMasterSlaveConstraints](../../functions/sources/communicator/Communicator/GlobalNumberOfMasterSlaveConstraints.md)
- [GetNumberOfColors](../../functions/sources/communicator/Communicator/GetNumberOfColors.md)
- [SetNumberOfColors](../../functions/sources/communicator/Communicator/SetNumberOfColors.md)
- [AddColors](../../functions/sources/communicator/Communicator/AddColors.md)
- [SetLocalMesh](../../functions/sources/communicator/Communicator/SetLocalMesh.md)
- [pLocalMesh](../../functions/sources/communicator/Communicator/pLocalMesh.md)
- [pGhostMesh](../../functions/sources/communicator/Communicator/pGhostMesh.md)
- [pInterfaceMesh](../../functions/sources/communicator/Communicator/pInterfaceMesh.md)
- [pLocalMesh](../../functions/sources/communicator/Communicator/pLocalMesh-2.md)
- [pGhostMesh](../../functions/sources/communicator/Communicator/pGhostMesh-2.md)
- [pInterfaceMesh](../../functions/sources/communicator/Communicator/pInterfaceMesh-2.md)
- [pLocalMesh](../../functions/sources/communicator/Communicator/pLocalMesh-3.md)
- [pGhostMesh](../../functions/sources/communicator/Communicator/pGhostMesh-3.md)
- [pInterfaceMesh](../../functions/sources/communicator/Communicator/pInterfaceMesh-3.md)
- [pLocalMesh](../../functions/sources/communicator/Communicator/pLocalMesh-4.md)
- [pGhostMesh](../../functions/sources/communicator/Communicator/pGhostMesh-4.md)
- [pInterfaceMesh](../../functions/sources/communicator/Communicator/pInterfaceMesh-4.md)
- [SynchronizeNodalSolutionStepsData](../../functions/sources/communicator/Communicator/SynchronizeNodalSolutionStepsData.md)
- [SynchronizeDofs](../../functions/sources/communicator/Communicator/SynchronizeDofs.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-2.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-3.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-4.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-5.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-6.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-7.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-8.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-9.md)
- [SynchronizeVariable](../../functions/sources/communicator/Communicator/SynchronizeVariable-10.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-2.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-3.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-4.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-5.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-6.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-7.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-8.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-9.md)
- [SynchronizeNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalVariable-10.md)
- [SynchronizeCurrentDataToMax](../../functions/sources/communicator/Communicator/SynchronizeCurrentDataToMax.md)
- [SynchronizeNonHistoricalDataToMax](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalDataToMax.md)
- [SynchronizeCurrentDataToAbsMax](../../functions/sources/communicator/Communicator/SynchronizeCurrentDataToAbsMax.md)
- [SynchronizeNonHistoricalDataToAbsMax](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalDataToAbsMax.md)
- [SynchronizeCurrentDataToMin](../../functions/sources/communicator/Communicator/SynchronizeCurrentDataToMin.md)
- [SynchronizeNonHistoricalDataToMin](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalDataToMin.md)
- [SynchronizeCurrentDataToAbsMin](../../functions/sources/communicator/Communicator/SynchronizeCurrentDataToAbsMin.md)
- [SynchronizeNonHistoricalDataToAbsMin](../../functions/sources/communicator/Communicator/SynchronizeNonHistoricalDataToAbsMin.md)
- [SynchronizeElementalFlags](../../functions/sources/communicator/Communicator/SynchronizeElementalFlags.md)
- [AssembleCurrentData](../../functions/sources/communicator/Communicator/AssembleCurrentData.md)
- [AssembleCurrentData](../../functions/sources/communicator/Communicator/AssembleCurrentData-2.md)
- [AssembleCurrentData](../../functions/sources/communicator/Communicator/AssembleCurrentData-3.md)
- [AssembleCurrentData](../../functions/sources/communicator/Communicator/AssembleCurrentData-4.md)
- [AssembleCurrentData](../../functions/sources/communicator/Communicator/AssembleCurrentData-5.md)
- [AssembleNonHistoricalData](../../functions/sources/communicator/Communicator/AssembleNonHistoricalData.md)
- [AssembleNonHistoricalData](../../functions/sources/communicator/Communicator/AssembleNonHistoricalData-2.md)
- [AssembleNonHistoricalData](../../functions/sources/communicator/Communicator/AssembleNonHistoricalData-3.md)
- [AssembleNonHistoricalData](../../functions/sources/communicator/Communicator/AssembleNonHistoricalData-4.md)
- [AssembleNonHistoricalData](../../functions/sources/communicator/Communicator/AssembleNonHistoricalData-5.md)
- [AssembleNonHistoricalData](../../functions/sources/communicator/Communicator/AssembleNonHistoricalData-6.md)
- [SynchronizeElementalNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeElementalNonHistoricalVariable.md)
- [SynchronizeElementalNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeElementalNonHistoricalVariable-2.md)
- [SynchronizeElementalNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeElementalNonHistoricalVariable-3.md)
- [SynchronizeElementalNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeElementalNonHistoricalVariable-4.md)
- [SynchronizeElementalNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeElementalNonHistoricalVariable-5.md)
- [SynchronizeElementalNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeElementalNonHistoricalVariable-6.md)
- [SynchronizeElementalNonHistoricalVariable](../../functions/sources/communicator/Communicator/SynchronizeElementalNonHistoricalVariable-7.md)
- [TransferObjects](../../functions/sources/communicator/Communicator/TransferObjects.md)
- [TransferObjects](../../functions/sources/communicator/Communicator/TransferObjects-2.md)
- [TransferObjects](../../functions/sources/communicator/Communicator/TransferObjects-3.md)
- [TransferObjects](../../functions/sources/communicator/Communicator/TransferObjects-4.md)
- [TransferObjects](../../functions/sources/communicator/Communicator/TransferObjects-5.md)
- [TransferObjects](../../functions/sources/communicator/Communicator/TransferObjects-6.md)
- [SynchronizeOrNodalFlags](../../functions/sources/communicator/Communicator/SynchronizeOrNodalFlags.md)
- [SynchronizeAndNodalFlags](../../functions/sources/communicator/Communicator/SynchronizeAndNodalFlags.md)
- [SynchronizeNodalFlags](../../functions/sources/communicator/Communicator/SynchronizeNodalFlags.md)
- [Clear](../../functions/sources/communicator/Communicator/Clear.md)
- [Info](../../functions/sources/communicator/Communicator/Info.md)
- [PrintInfo](../../functions/sources/communicator/Communicator/PrintInfo.md)
- [PrintData](../../functions/sources/communicator/Communicator/PrintData.md)

# Imports

- `"includes/communicator.h"`
- `"includes/data_communicator.h"`
- `"includes/parallel_environment.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)