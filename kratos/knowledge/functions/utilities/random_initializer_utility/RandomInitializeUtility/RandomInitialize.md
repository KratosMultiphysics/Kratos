---
type: C++ Method
title: RandomInitialize
resource: utilities/random_initializer_utility.h#L118-L128
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/random_initializer_utility/RandomInitializeUtility/NormalDestributionRandom
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/linear_solvers/power_iteration_eigenvalue_solver/PowerIterationEigenvalueSolver/Solve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/power_iteration_highest_eigenvalue_solver/PowerIterationHighestEigenvalueSolver/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void RandomInitialize( const SparseMatrixType& K, VectorType& R, const bool Inverse = false )`

# Calls

- [NormalDestributionRandom](../../../../functions/utilities/random_initializer_utility/RandomInitializeUtility/NormalDestributionRandom.md)

# Called by

- [Solve](../../../../functions/linear_solvers/power_iteration_eigenvalue_solver/PowerIterationEigenvalueSolver/Solve.md)
- [Solve](../../../../functions/linear_solvers/power_iteration_highest_eigenvalue_solver/PowerIterationHighestEigenvalueSolver/Solve.md)