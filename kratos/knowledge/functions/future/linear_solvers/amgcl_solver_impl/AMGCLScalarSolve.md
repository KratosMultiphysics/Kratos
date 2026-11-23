---
type: C++ Function
title: AMGCLScalarSolve
resource: future/linear_solvers/amgcl_solver_impl.cpp#L43-L104
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/linear_solvers/amgcl_solver_impl/AMGCLSolve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AMGCLScalarSolve( CsrMatrix<>& rA, SystemVector<>& rX, SystemVector<>& rB, typename CsrMatrix<>::IndexType& rIterationNumber, double& rResidual, const boost::property_tree::ptree &amgclParams, int verbosity_level, bool use_gpgpu )`

# Called by

- [AMGCLSolve](../../../../functions/future/linear_solvers/amgcl_solver_impl/AMGCLSolve.md)