---
type: C++ Method
title: AMGCLSolver
resource: linear_solvers/amgcl_solver_impl.hpp#L189-L213
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`AMGCLSolver<TSparse,TDense>::AMGCLSolver(const std::string& rSmootherName, const std::string& rSolverName, double Tolerance, int MaxIterationsNumber, int Verbosity, int GMRESSize) : mTolerance(Tolerance), mMaxIterationsNumber(MaxIterationsNumber), mVerbosity(Verbosity), mBlockSize(), mGMRESSize(GMRESSize), mCoarseEnough(1000), mFallbackToGMRES(false), mProvideCoordinates(false), mUseBlockMatricesIfPossible(false)`