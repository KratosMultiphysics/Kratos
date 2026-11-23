---
type: C++ Module
title: simple_mortar_mapper_process
resource: processes/simple_mortar_mapper_process.cpp#L1-L1225
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/processes-simple-mortar-mapper-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometrical-projection-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-mortar-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-normal-calculation-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-math-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-atomic-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [SimpleMortarMapperProcess](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/SimpleMortarMapperProcess.md)
- [SimpleMortarMapperProcess](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/SimpleMortarMapperProcess-2.md)
- [SimpleMortarMapperProcess](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/SimpleMortarMapperProcess-3.md)
- [Execute](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/Execute.md)
- [ExecuteInitializeSolutionStep](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteInitializeSolutionStep.md)
- [CheckAndPerformSearch](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CheckAndPerformSearch.md)
- [ResetNodalArea](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ResetNodalArea.md)
- [GetReferenceArea](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetReferenceArea.md)
- [AssemblyMortarOperators](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators.md)
- [CalculateAe](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CalculateAe.md)
- [InvertDiagonalMatrix](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/InvertDiagonalMatrix.md)
- [InvertDiagonalMatrix](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/InvertDiagonalMatrix-2.md)
- [LumpMatrix](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/LumpMatrix.md)
- [GetSystemSize](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetSystemSize.md)
- [CreateSlaveConectivityDatabase](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CreateSlaveConectivityDatabase.md)
- [GetIntegrationMethod](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetIntegrationMethod.md)
- [CheckWholeVector](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CheckWholeVector.md)
- [ComputeResidualMatrix](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ComputeResidualMatrix.md)
- [AssembleRHSAndLHS](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssembleRHSAndLHS.md)
- [AssembleRHS](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssembleRHS.md)
- [ExecuteExplicitMapping](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteExplicitMapping.md)
- [ExecuteImplicitMapping](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping.md)
- [CreateInverseDatabase](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CreateInverseDatabase.md)
- [UpdateInterface](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/UpdateInterface.md)
- [GetDefaultParameters](../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetDefaultParameters.md)

# Imports

- `"processes/simple_mortar_mapper_process.h"`
- `"utilities/geometrical_projection_utilities.h"`
- `"utilities/mortar_utilities.h"`
- `"utilities/normal_calculation_utils.h"`
- `"utilities/math_utils.h"`
- `"utilities/atomic_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)