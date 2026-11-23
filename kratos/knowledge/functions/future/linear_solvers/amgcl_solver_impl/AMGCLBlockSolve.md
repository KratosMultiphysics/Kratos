---
type: C++ Function
title: AMGCLBlockSolve
resource: future/linear_solvers/amgcl_solver_impl.cpp#L107-L187
generated:
  by: okf-rs/0.7.0
---

# Signature

`void AMGCLBlockSolve( CsrMatrix<> & rA, SystemVector<>& rX, SystemVector<>& rB, typename CsrMatrix<>::IndexType& rIterationNumber, double& rResidual, boost::property_tree::ptree amgclParams, int verbosity_level, bool use_gpgpu )`