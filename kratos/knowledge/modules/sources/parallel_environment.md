---
type: C++ Module
title: parallel_environment
resource: sources/parallel_environment.cpp#L1-L396
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/mutex
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-parallel-environment-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-components-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/input-output-logger-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [SetDefaultDataCommunicator](../../functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicator.md)
- [GetDefaultRank](../../functions/sources/parallel_environment/ParallelEnvironment/GetDefaultRank.md)
- [GetDefaultSize](../../functions/sources/parallel_environment/ParallelEnvironment/GetDefaultSize.md)
- [SetUpMPIEnvironment](../../functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironment.md)
- [RegisterFillCommunicatorFactory](../../functions/sources/parallel_environment/ParallelEnvironment/RegisterFillCommunicatorFactory.md)
- [RegisterCommunicatorFactory](../../functions/sources/parallel_environment/ParallelEnvironment/RegisterCommunicatorFactory.md)
- [CreateFillCommunicator](../../functions/sources/parallel_environment/ParallelEnvironment/CreateFillCommunicator.md)
- [CreateFillCommunicatorFromGlobalParallelism](../../functions/sources/parallel_environment/ParallelEnvironment/CreateFillCommunicatorFromGlobalParallelism.md)
- [CreateFillCommunicatorFromGlobalParallelism](../../functions/sources/parallel_environment/ParallelEnvironment/CreateFillCommunicatorFromGlobalParallelism-2.md)
- [CreateCommunicatorFromGlobalParallelism](../../functions/sources/parallel_environment/ParallelEnvironment/CreateCommunicatorFromGlobalParallelism.md)
- [CreateCommunicatorFromGlobalParallelism](../../functions/sources/parallel_environment/ParallelEnvironment/CreateCommunicatorFromGlobalParallelism-2.md)
- [RegisterDataCommunicator](../../functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicator.md)
- [UnregisterDataCommunicator](../../functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator.md)
- [HasDataCommunicator](../../functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicator.md)
- [GetDefaultDataCommunicatorName](../../functions/sources/parallel_environment/ParallelEnvironment/GetDefaultDataCommunicatorName.md)
- [MPIIsInitialized](../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsInitialized.md)
- [MPIIsFinalized](../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsFinalized.md)
- [Info](../../functions/sources/parallel_environment/ParallelEnvironment/Info.md)
- [PrintInfo](../../functions/sources/parallel_environment/ParallelEnvironment/PrintInfo.md)
- [PrintData](../../functions/sources/parallel_environment/ParallelEnvironment/PrintData.md)
- [RegisterFillCommunicatorFactoryDetail](../../functions/sources/parallel_environment/ParallelEnvironment/RegisterFillCommunicatorFactoryDetail.md)
- [RegisterFillCommunicatorFactoryDetail](../../functions/sources/parallel_environment/ParallelEnvironment/RegisterFillCommunicatorFactoryDetail-2.md)
- [RegisterCommunicatorFactoryDetail](../../functions/sources/parallel_environment/ParallelEnvironment/RegisterCommunicatorFactoryDetail.md)
- [RegisterCommunicatorFactoryDetail](../../functions/sources/parallel_environment/ParallelEnvironment/RegisterCommunicatorFactoryDetail-2.md)
- [ParallelEnvironment](../../functions/sources/parallel_environment/ParallelEnvironment/ParallelEnvironment.md)
- [RetrieveRegisteredName](../../functions/sources/parallel_environment/ParallelEnvironment/RetrieveRegisteredName.md)
- [SetAsDefault](../../functions/sources/parallel_environment/ParallelEnvironment/SetAsDefault.md)
- [Create](../../functions/sources/parallel_environment/ParallelEnvironment/Create.md)
- [SetUpMPIEnvironmentDetail](../../functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironmentDetail.md)
- [RegisterDataCommunicatorDetail](../../functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicatorDetail.md)
- [UnregisterDataCommunicatorDetail](../../functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicatorDetail.md)
- [SetDefaultDataCommunicatorDetail](../../functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicatorDetail.md)
- [HasDataCommunicatorDetail](../../functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicatorDetail.md)
- [MPIIsInitializedDetail](../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsInitializedDetail.md)
- [MPIIsFinalizedDetail](../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsFinalizedDetail.md)
- [InfoDetail](../../functions/sources/parallel_environment/ParallelEnvironment/InfoDetail.md)
- [PrintInfoDetail](../../functions/sources/parallel_environment/ParallelEnvironment/PrintInfoDetail.md)
- [PrintDataDetail](../../functions/sources/parallel_environment/ParallelEnvironment/PrintDataDetail.md)

# Imports

- `<mutex>`
- `"includes/parallel_environment.h"`
- `"includes/kratos_components.h"`
- `"input_output/logger.h"`
- `"utilities/parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)