---
type: C++ Module
title: residualbased_block_builder_and_solver
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L1-L1789
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/unordered-set
    resolved_by: tree-sitter
    confidence: exact
  - target: external/omp-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-builder-and-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-key-hash-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-lock-object-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-dof-utilities-block-build-dof-array-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-sparse-matrix-multiplication-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-builtin-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-atomic-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ResidualBasedBlockBuilderAndSolver](../../../classes/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver.md)
- [ResidualBasedBlockBuilderAndSolver](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ResidualBasedBlockBuilderAndSolver.md)
- [ResidualBasedBlockBuilderAndSolver](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ResidualBasedBlockBuilderAndSolver-2.md)
- [ResidualBasedBlockBuilderAndSolver](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ResidualBasedBlockBuilderAndSolver-3.md)
- [Create](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Create.md)
- [Build](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Build.md)
- [BuildLHS](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildLHS.md)
- [BuildLHS_CompleteOnFreeRows](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildLHS_CompleteOnFreeRows.md)
- [SystemSolve](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SystemSolve.md)
- [SystemSolveWithPhysics](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SystemSolveWithPhysics.md)
- [InternalSystemSolveWithPhysics](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/InternalSystemSolveWithPhysics.md)
- [BuildAndSolve](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolve.md)
- [BuildAndSolveLinearizedOnPreviousIteration](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration.md)
- [BuildRHSAndSolve](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHSAndSolve.md)
- [BuildRHS](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHS.md)
- [SetUpDofSet](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SetUpDofSet.md)
- [SetUpSystem](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SetUpSystem.md)
- [ResizeAndInitializeVectors](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ResizeAndInitializeVectors.md)
- [CalculateReactions](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/CalculateReactions.md)
- [ApplyDirichletConditions](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyDirichletConditions.md)
- [ApplyRHSConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyRHSConstraints.md)
- [ApplyConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyConstraints.md)
- [Clear](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Clear.md)
- [Check](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Check.md)
- [GetDefaultParameters](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/GetDefaultParameters.md)
- [Name](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Name.md)
- [GetScaleFactor](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/GetScaleFactor.md)
- [SetScaleFactor](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SetScaleFactor.md)
- [IsConstantConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/IsConstantConstraints.md)
- [SetConstantConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SetConstantConstraints.md)
- [Info](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Info.md)
- [PrintInfo](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/PrintInfo.md)
- [PrintData](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/PrintData.md)
- [BuildRHSNoDirichlet](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHSNoDirichlet.md)
- [ConstructMasterSlaveConstraintsStructure](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ConstructMasterSlaveConstraintsStructure.md)
- [BuildMasterSlaveConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildMasterSlaveConstraints.md)
- [for](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/for.md)
- [ConstructMatrixStructure](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ConstructMatrixStructure.md)
- [TLS](../../../classes/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/TLS.md)
- [Assemble](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Assemble.md)
- [AssembleLHS](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/AssembleLHS.md)
- [AssembleRHS](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/AssembleRHS.md)
- [AssembleRowContribution](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/AssembleRowContribution.md)
- [AssignSettings](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/AssignSettings.md)
- [AddUnique](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/AddUnique.md)
- [CreatePartition](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/CreatePartition.md)
- [ForwardFind](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ForwardFind.md)
- [BackwardFind](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BackwardFind.md)
- [GetCurrentThreadId](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/GetCurrentThreadId.md)

# Imports

- `<unordered_set>`
- `<omp.h>`
- `"solving_strategies/builder_and_solvers/builder_and_solver.h"`
- `"includes/model_part.h"`
- `"includes/key_hash.h"`
- `"utilities/timer.h"`
- `"utilities/variable_utils.h"`
- `"includes/kratos_flags.h"`
- `"includes/lock_object.h"`
- `"utilities/dof_utilities/block_build_dof_array_utility.h"`
- `"utilities/sparse_matrix_multiplication_utility.h"`
- `"utilities/builtin_timer.h"`
- `"utilities/atomic_utilities.h"`
- `"spaces/ublas_space.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)