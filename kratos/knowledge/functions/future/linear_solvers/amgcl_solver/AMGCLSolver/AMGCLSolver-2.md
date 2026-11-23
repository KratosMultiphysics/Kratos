---
type: C++ Method
title: AMGCLSolver
resource: future/linear_solvers/amgcl_solver.h#L246-L271
generated:
  by: okf-rs/0.7.0
---

# Signature

`AMGCLSolver( AMGCLSmoother Smoother, AMGCLIterativeSolverType Solver, double Tolerance, int MaxIterationsNumber, int Verbosity, int GMRESSize = 50 ) : mTolerance(Tolerance), mMaxIterationsNumber(MaxIterationsNumber), mVerbosity(Verbosity), mBlockSize(1), mGMRESSize(GMRESSize), mCoarseEnough(1000), mFallbackToGMRES(false), mProvideCoordinates(false), mUseBlockMatricesIfPossible(false)`