---
type: C++ Module
title: augmented_lagrange_constraint_assembler
resource: solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.cpp#L1-L517
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-augmented-lagrange-constraint-assembler-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-sparse-utilities-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-constraint-utilities-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-diagonal-scaling-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-sparse-matrix-multiplication-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/filesystem
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ApplyDirichletConditions](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/ApplyDirichletConditions.md)
- [AugmentedLagrangeConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AugmentedLagrangeConstraintAssembler.md)
- [AugmentedLagrangeConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AugmentedLagrangeConstraintAssembler-2.md)
- [AugmentedLagrangeConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AugmentedLagrangeConstraintAssembler-3.md)
- [AllocateConstraints](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AllocateConstraints.md)
- [AllocateSystem](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AllocateSystem.md)
- [Assemble](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Assemble.md)
- [Initialize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Initialize.md)
- [InitializeConstraintIteration](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/InitializeConstraintIteration.md)
- [FinalizeConstraintIteration](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/FinalizeConstraintIteration.md)
- [Finalize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Finalize.md)
- [Clear](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Clear.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/GetDefaultParameters.md)
- [GetPenaltyFactor](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/GetPenaltyFactor.md)

# Imports

- `"solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp"`
- `"spaces/ublas_space.h"`
- `"utilities/sparse_matrix_multiplication_utility.h"`
- `<filesystem>`

# Member of

- [KratosMultiphysics](../../../../packages/KratosMultiphysics.md)