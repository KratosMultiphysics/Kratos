---
type: C++ Method
title: AMGCLSolver
resource: linear_solvers/amgcl_solver_impl.hpp#L217-L243
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`AMGCLSolver<TSparse,TDense>::AMGCLSolver(const std::string& rSmootherName, const std::string& rSolverName, const std::string& rCoarseningName, double Tolerance, int MaxIterationsNumber, int Verbosity, int GMRESSize, bool ProvideCoordinates) : mTolerance(Tolerance), mMaxIterationsNumber(MaxIterationsNumber), mVerbosity(Verbosity), mBlockSize(), mGMRESSize(GMRESSize), mCoarseEnough(1000), mFallbackToGMRES(false), mProvideCoordinates(ProvideCoordinates), mUseBlockMatricesIfPossible(false)`