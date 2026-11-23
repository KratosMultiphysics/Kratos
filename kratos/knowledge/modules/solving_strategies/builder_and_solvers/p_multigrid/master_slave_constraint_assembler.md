---
type: C++ Module
title: master_slave_constraint_assembler
resource: solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler.cpp#L1-L628
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-master-slave-constraint-assembler-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-diagonal-scaling-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-sparse-utilities-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-constraint-utilities-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/factories-linear-solver-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-sparse-matrix-multiplication-utility-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [LinearSystem](../../../../classes/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/LinearSystem.md)
- [GetIndependentTransform](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Impl/GetIndependentTransform.md)
- [GetDependentTransform](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Impl/GetDependentTransform.md)
- [MasterSlaveConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/MasterSlaveConstraintAssembler.md)
- [MasterSlaveConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/MasterSlaveConstraintAssembler-2.md)
- [MasterSlaveConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/MasterSlaveConstraintAssembler-3.md)
- [AllocateConstraints](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/AllocateConstraints.md)
- [Assemble](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Assemble.md)
- [Initialize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Initialize.md)
- [FinalizeConstraintIteration](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/FinalizeConstraintIteration.md)
- [Finalize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Finalize.md)
- [ComputeIndependentResidual](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/ComputeIndependentResidual.md)
- [ComputeDependentResidual](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/ComputeDependentResidual.md)
- [ComputeIndependentSolution](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/ComputeIndependentSolution.md)
- [ComputeDependentSolution](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/ComputeDependentSolution.md)
- [Clear](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Clear.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/GetDefaultParameters.md)

# Imports

- `"solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities.hpp"`
- `"factories/linear_solver_factory.h"`
- `"includes/define.h"`
- `"spaces/ublas_space.h"`
- `"utilities/sparse_matrix_multiplication_utility.h"`

# Member of

- [KratosMultiphysics](../../../../packages/KratosMultiphysics.md)