---
type: C++ Method
title: IterativeSolve
resource: linear_solvers/deflated_cg_solver.h#L320-L453
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/deflation_utils/FillDeflatedMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/skyline_lu_factorization_solver/LUSkylineFactorization/copyFromCSRMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/deflation_utils/ApplyWtranspose
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/deflation_utils/ApplyW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/iterative_solver/IterativeSolver/PreconditionedMult
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IterativeSolve(SparseMatrixType& rA, SparseVectorType& rX, SparseVectorType& rB)`

# Calls

- [FillDeflatedMatrix](../../../../functions/utilities/deflation_utils/FillDeflatedMatrix.md)
- [copyFromCSRMatrix](../../../../functions/future/linear_solvers/skyline_lu_factorization_solver/LUSkylineFactorization/copyFromCSRMatrix.md)
- [ApplyWtranspose](../../../../functions/utilities/deflation_utils/ApplyWtranspose.md)
- [ApplyW](../../../../functions/utilities/deflation_utils/ApplyW.md)
- [PreconditionedMult](../../../../functions/linear_solvers/iterative_solver/IterativeSolver/PreconditionedMult.md)