---
type: C++ Method
title: nnz
resource: containers/csr_matrix.h#L269-L272
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/csr_matrix/SpMV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/ApplyHomogeneousDirichlet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/local_nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ApplyHomogeneousDirichlet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToSerialCSR
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/AllocateConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/InPlaceMatrixAdd
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/MergeMatrices
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/CheckMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/BalancedProduct
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ResetConstraintSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/matrix_market_interface/WriteMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/sparse_containers_test_utilities/SparseContainersTestUtilities/CheckCSRMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/deflation_utils/ConstructW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/AssembleSparseMatrixByBlocks
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline IndexType nnz() const`

# Called by

- [SpMV](../../../../functions/containers/csr_matrix/SpMV.md)
- [ApplyHomogeneousDirichlet](../../../../functions/containers/csr_matrix/ApplyHomogeneousDirichlet.md)
- [PrintData](../../../../functions/containers/csr_matrix/CsrMatrix/PrintData.md)
- [local_nnz](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/local_nnz.md)
- [ApplyHomogeneousDirichlet](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ApplyHomogeneousDirichlet.md)
- [ToSerialCSR](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToSerialCSR.md)
- [PrintData](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrintData.md)
- [Initialize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Initialize.md)
- [AllocateConstraints](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/AllocateConstraints.md)
- [InPlaceMatrixAdd](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/InPlaceMatrixAdd.md)
- [MergeMatrices](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/MergeMatrices.md)
- [CheckMatrix](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/CheckMatrix.md)
- [BalancedProduct](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/BalancedProduct.md)
- [ApplyConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints.md)
- [ResetConstraintSystem](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ResetConstraintSystem.md)
- [WriteMatrixMarketMatrix](../../../../functions/sources/matrix_market_interface/WriteMatrixMarketMatrix.md)
- [CheckCSRMatrix](../../../../functions/tests/test_utilities/sparse_containers_test_utilities/SparseContainersTestUtilities/CheckCSRMatrix.md)
- [ConstructW](../../../../functions/utilities/deflation_utils/ConstructW.md)
- [AssembleSparseMatrixByBlocks](../../../../functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/AssembleSparseMatrixByBlocks.md)