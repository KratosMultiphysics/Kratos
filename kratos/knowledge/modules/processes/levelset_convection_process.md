---
type: C++ Module
title: levelset_convection_process
resource: processes/levelset_convection_process.h#L1-L1124
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/algorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-convection-diffusion-settings-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-global-pointer-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/elements-levelset-convection-element-simplex-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/elements-levelset-convection-element-simplex-algebraic-stabilization-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-geometry-data-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-compute-nodal-gradient-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-schemes-residualbased-incrementalupdate-static-scheme-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-residualbased-block-builder-and-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-strategies-residualbased-linear-strategy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-pointer-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-pointer-map-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-nodal-h-process-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [LevelSetConvectionProcess](../../classes/processes/levelset_convection_process/LevelSetConvectionProcess.md)
- [ProcessInfoDataContainer](../../classes/processes/levelset_convection_process/ProcessInfoDataContainer.md)
- [ProcessInfoDataContainer](../../functions/processes/levelset_convection_process/ProcessInfoDataContainer/ProcessInfoDataContainer.md)
- [RestoreProcessInfoData](../../functions/processes/levelset_convection_process/ProcessInfoDataContainer/RestoreProcessInfoData.md)
- [LevelSetConvectionProcess](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess.md)
- [LevelSetConvectionProcess](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess-2.md)
- [LevelSetConvectionProcess](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess-3.md)
- [Execute](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute.md)
- [Clear](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Clear.md)
- [GetDefaultParameters](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetDefaultParameters.md)
- [Info](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Info.md)
- [PrintInfo](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/PrintInfo.md)
- [PrintData](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/PrintData.md)
- [LevelSetConvectionProcess](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess-4.md)
- [SetConvectionProblemSettings](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/SetConvectionProblemSettings.md)
- [ReGenerateConvectionModelPart](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/ReGenerateConvectionModelPart.md)
- [InitializeDistanceModelPartDatabases](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeDistanceModelPartDatabases.md)
- [EvaluateNumberOfSubsteps](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/EvaluateNumberOfSubsteps.md)
- [EvaluateLimiter](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/EvaluateLimiter.md)
- [ErrorCalculationAndCorrection](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/ErrorCalculationAndCorrection.md)
- [ComputeNodalH](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/ComputeNodalH.md)
- [CheckAndAssignSettings](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/CheckAndAssignSettings.md)
- [GetConvectionElementRegisteredName](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementRegisteredName.md)
- [GetConvectionElementsList](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementsList.md)
- [GetConvectionElementName](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementName.md)
- [GetConvectionElementDefaultParameters](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementDefaultParameters.md)
- [GetFillProcessInfoFormulationDataFunction](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetFillProcessInfoFormulationDataFunction.md)
- [InitializeConvectionStrategy](../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeConvectionStrategy.md)

# Imports

- `<string>`
- `<iostream>`
- `<algorithm>`
- `"includes/convection_diffusion_settings.h"`
- `"includes/define.h"`
- `"includes/global_pointer_variables.h"`
- `"includes/kratos_flags.h"`
- `"elements/levelset_convection_element_simplex.h"`
- `"elements/levelset_convection_element_simplex_algebraic_stabilization.h"`
- `"geometries/geometry_data.h"`
- `"processes/compute_nodal_gradient_process.h"`
- `"solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"`
- `"solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"`
- `"solving_strategies/strategies/residualbased_linear_strategy.h"`
- `"utilities/variable_utils.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/pointer_communicator.h"`
- `"utilities/pointer_map_communicator.h"`
- `"processes/find_nodal_h_process.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)