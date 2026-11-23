---
type: C++ Function
title: AMGCLSolve
resource: future/linear_solvers/amgcl_solver_impl.cpp#L189-L225
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/linear_solvers/amgcl_solver_impl/AMGCLScalarSolve
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/linear_solvers/amgcl_solver/AMGCLSolver/PerformSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AMGCLSolve( int block_size, CsrMatrix<>& rA, SystemVector<>& rX, SystemVector<>& rB, typename CsrMatrix<>::IndexType& rIterationNumber, double& rResidual, boost::property_tree::ptree amgclParams, int verbosity_level, bool use_gpgpu )`

# Calls

- [AMGCLScalarSolve](../../../../functions/future/linear_solvers/amgcl_solver_impl/AMGCLScalarSolve.md)

# Called by

- [PerformSolutionStep](../../../../functions/future/linear_solvers/amgcl_solver/AMGCLSolver/PerformSolutionStep.md)